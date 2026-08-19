#!/usr/bin/env python3
"""Part 4 -- baseline subset analysis.

Subset of interest: melanoma patients, PBMC samples, treated with miraclib,
at baseline (time_from_treatment_start = 0).  For that subset we report

  a) samples per project,
  b) responders vs non-responders,
  c) males vs females,

and finally the average B-cell count for melanoma **males who responded** at
time = 0 across *all* sample types and treatments.

Everything here is answered in SQL against the normalised schema -- the point
of Part 1 is that these questions become filters, not bespoke code.
"""
from __future__ import annotations

import argparse

import pandas as pd

from analysis.common import banner, query, write_csv

BASELINE_SUBSET = """
SELECT  s.sample_code   AS sample,
        sub.subject_code AS subject,
        p.project_code   AS project,
        e.response       AS response,
        sub.sex          AS sex,
        e.age_at_enrollment AS age
FROM        sample      s
JOIN        enrollment  e   ON e.enrollment_id   = s.enrollment_id
JOIN        subject     sub ON sub.subject_id    = e.subject_id
JOIN        project     p   ON p.project_id      = e.project_id
JOIN        condition   c   ON c.condition_id    = e.condition_id
JOIN        treatment   t   ON t.treatment_id    = e.treatment_id
JOIN        sample_type st  ON st.sample_type_id = s.sample_type_id
WHERE   c.condition_name   = 'melanoma'
  AND   t.treatment_name   = 'miraclib'
  AND   st.sample_type_name = 'PBMC'
  AND   s.time_from_treatment_start = 0
"""

SAMPLES_PER_PROJECT = """
SELECT project, COUNT(*) AS n_samples,
       COUNT(DISTINCT subject) AS n_subjects
FROM   ({base})
GROUP  BY project
ORDER  BY project
""".format(base=BASELINE_SUBSET)

RESPONSE_BREAKDOWN = """
SELECT COALESCE(response, 'unknown') AS response,
       COUNT(*) AS n_samples,
       COUNT(DISTINCT subject) AS n_subjects
FROM   ({base})
GROUP  BY COALESCE(response, 'unknown')
ORDER  BY response DESC
""".format(base=BASELINE_SUBSET)

SEX_BREAKDOWN = """
SELECT COALESCE(sex, 'unknown') AS sex,
       COUNT(*) AS n_samples,
       COUNT(DISTINCT subject) AS n_subjects
FROM   ({base})
GROUP  BY COALESCE(sex, 'unknown')
ORDER  BY sex
""".format(base=BASELINE_SUBSET)

# All sample types, all treatments -- melanoma males, responders, baseline.
AVG_B_CELL = """
SELECT  COUNT(*)            AS n_samples,
        COUNT(DISTINCT sub.subject_code) AS n_subjects,
        AVG(scc.cell_count) AS avg_b_cell_count,
        MIN(scc.cell_count) AS min_b_cell_count,
        MAX(scc.cell_count) AS max_b_cell_count
FROM        sample_cell_count scc
JOIN        cell_population cp  ON cp.population_id = scc.population_id
JOIN        sample      s   ON s.sample_id       = scc.sample_id
JOIN        enrollment  e   ON e.enrollment_id   = s.enrollment_id
JOIN        subject     sub ON sub.subject_id    = e.subject_id
JOIN        condition   c   ON c.condition_id    = e.condition_id
WHERE   cp.population_name = 'b_cell'
  AND   c.condition_name   = 'melanoma'
  AND   sub.sex            = 'M'
  AND   e.response         = 'yes'
  AND   s.time_from_treatment_start = 0
"""


def baseline_subset(db_path: str = None) -> pd.DataFrame:
    return query(BASELINE_SUBSET + " ORDER BY sample", db_path=db_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db", default=None)
    args = parser.parse_args()

    subset = baseline_subset(args.db)
    by_project = query(SAMPLES_PER_PROJECT, db_path=args.db)
    by_response = query(RESPONSE_BREAKDOWN, db_path=args.db)
    by_sex = query(SEX_BREAKDOWN, db_path=args.db)
    b_cell = query(AVG_B_CELL, db_path=args.db)

    banner("PART 4 -- Melanoma / miraclib / PBMC / baseline (time = 0)")
    print("Subset size: %d samples from %d subjects\n"
          % (len(subset), subset["subject"].nunique()))

    print("(a) Samples per project")
    print(by_project.to_string(index=False))
    print("\n(b) Responders vs non-responders")
    print(by_response.to_string(index=False))
    print("\n(c) Males vs females")
    print(by_sex.to_string(index=False))

    banner("PART 4 -- Average B cells, melanoma males, responders, time = 0"
           " (all sample types & treatments)")
    row = b_cell.iloc[0]
    if pd.isna(row["avg_b_cell_count"]):
        print("No matching samples.")
    else:
        print("n = %d samples (%d subjects)" % (row["n_samples"], row["n_subjects"]))
        print("average B cell count = %.2f  (range %d - %d)"
              % (row["avg_b_cell_count"], row["min_b_cell_count"], row["max_b_cell_count"]))

    write_csv(subset, "part4_baseline_subset.csv")
    write_csv(by_project, "part4_samples_per_project.csv")
    write_csv(by_response, "part4_response_breakdown.csv")
    write_csv(by_sex, "part4_sex_breakdown.csv")
    write_csv(b_cell, "part4_avg_b_cell_male_responders.csv")
    print("\nTables written to outputs/part4_*.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

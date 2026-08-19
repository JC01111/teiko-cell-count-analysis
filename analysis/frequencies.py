#!/usr/bin/env python3
"""Part 2 -- relative frequency of each cell population per sample.

Produces a long-format table with one row per (sample, population):

    sample | total_count | population | count | percentage

`total_count` is the sum of all five populations for that sample and
`percentage` is 100 * count / total_count.  The arithmetic lives in the
`v_sample_frequency` SQL view so the dashboard, the stats module and this
script can never disagree about how a frequency is defined.
"""
from __future__ import annotations

import argparse

import pandas as pd

from analysis.common import banner, query, write_csv

SQL = """
SELECT sample, total_count, population, count, percentage
FROM   v_sample_frequency
ORDER  BY sample, population
"""


def cell_frequencies(db_path: str = None) -> pd.DataFrame:
    return query(SQL, db_path=db_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db", default=None)
    parser.add_argument("--out", default="cell_frequencies.csv")
    args = parser.parse_args()

    df = cell_frequencies(args.db)
    path = write_csv(df, args.out)

    banner("PART 2 -- Relative frequency of each cell population per sample")
    print(df.head(10).to_string(index=False))
    print("...")
    print("\n%d rows (%d samples x %d populations) -> %s"
          % (len(df), df["sample"].nunique(), df["population"].nunique(), path))

    # Sanity check: percentages must sum to 100 within every sample.
    sums = df.groupby("sample")["percentage"].sum()
    print("Per-sample percentage sum: min=%.4f max=%.4f" % (sums.min(), sums.max()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Smoke + correctness tests for the loader and the analysis modules."""
from __future__ import annotations

import os
import sqlite3
import textwrap

import pytest

import load_data
from analysis import frequencies, statistics_report, subsets
from analysis.common import analysis_frame

CSV = textwrap.dedent(
    """\
    project,subject,condition,age,sex,treatment,response,sample,sample_type,time_from_treatment_start,b_cell,cd8_t_cell,cd4_t_cell,nk_cell,monocyte
    prj1,sbj000,melanoma,57,M,miraclib,yes,s0,PBMC,0,100,200,300,200,200
    prj1,sbj000,melanoma,57,M,miraclib,yes,s1,PBMC,7,150,150,300,200,200
    prj1,sbj001,melanoma,68,F,miraclib,no,s2,PBMC,0,300,200,200,150,150
    prj2,sbj002,healthy,44,F,none,,s3,WB,0,100,100,100,100,100
    """
)


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    csv_path = d / "cell-count.csv"
    csv_path.write_text(CSV)
    db_path = d / "test.db"
    assert load_data.main(["--csv", str(csv_path), "--db", str(db_path)]) == 0
    return str(db_path)


def test_row_counts(db):
    conn = sqlite3.connect(db)
    counts = {t: conn.execute("SELECT COUNT(*) FROM %s" % t).fetchone()[0]
              for t in ("subject", "enrollment", "sample", "sample_cell_count",
                        "cell_population", "project", "condition", "treatment", "sample_type")}
    conn.close()
    assert counts["subject"] == 3
    assert counts["enrollment"] == 3
    assert counts["sample"] == 4
    assert counts["sample_cell_count"] == 20
    assert counts["cell_population"] == 5
    assert counts["treatment"] == 2  # 'none' is a real treatment, not a NULL


def test_healthy_response_is_null(db):
    conn = sqlite3.connect(db)
    val = conn.execute(
        "SELECT e.response FROM enrollment e JOIN subject s USING (subject_id)"
        " WHERE s.subject_code = 'sbj002'").fetchone()[0]
    conn.close()
    assert val is None


def test_frequencies_sum_to_100(db):
    df = frequencies.cell_frequencies(db)
    assert len(df) == 20
    sums = df.groupby("sample")["percentage"].sum()
    assert sums.between(99.99, 100.01).all()
    row = df[(df["sample"] == "s0") & (df["population"] == "b_cell")].iloc[0]
    assert row["total_count"] == 1000
    assert row["percentage"] == pytest.approx(10.0)


def test_reload_is_idempotent(db, tmp_path):
    csv_path = tmp_path / "c.csv"
    csv_path.write_text(CSV)
    before = sqlite3.connect(db).execute("SELECT COUNT(*) FROM sample_cell_count").fetchone()[0]
    load_data.main(["--csv", str(csv_path), "--db", db, "--append"])
    after = sqlite3.connect(db).execute("SELECT COUNT(*) FROM sample_cell_count").fetchone()[0]
    assert before == after


def test_part3_cohort_filter(db):
    df = analysis_frame(db)
    cohort = statistics_report.melanoma_miraclib_pbmc(df)
    assert set(cohort["sample"]) == {"s0", "s1", "s2"}
    assert set(cohort["condition"]) == {"melanoma"}


def test_bh_fdr_monotonic():
    import numpy as np
    p = np.array([0.001, 0.01, 0.04, 0.5])
    adj = statistics_report._bh_fdr(p)
    assert (adj >= p).all()
    assert (np.diff(adj) >= -1e-12).all()


def test_part4_subset(db):
    sub = subsets.baseline_subset(db)
    assert set(sub["sample"]) == {"s0", "s2"}   # PBMC + melanoma + miraclib + t=0
    assert sub["project"].tolist() == ["prj1", "prj1"]


def test_new_population_column_needs_no_schema_change(tmp_path):
    csv_path = tmp_path / "extra.csv"
    csv_path.write_text(CSV.replace("monocyte\n", "monocyte,dendritic_cell\n")
                        .replace(",200,200\n", ",200,200,50\n")
                        .replace(",150,150\n", ",150,150,50\n")
                        .replace(",100,100\n", ",100,100,50\n"))
    db_path = tmp_path / "extra.db"
    assert load_data.main(["--csv", str(csv_path), "--db", str(db_path)]) == 0
    conn = sqlite3.connect(str(db_path))
    pops = [r[0] for r in conn.execute("SELECT population_name FROM cell_population")]
    conn.close()
    assert "dendritic_cell" in pops

"""Shared paths, database access and small helpers for the analysis scripts."""
from __future__ import annotations

import os
import sqlite3
from typing import Optional

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.environ.get("CELLCOUNT_DB", os.path.join(ROOT, "db", "cell_counts.db"))
OUTPUT_DIR = os.path.join(ROOT, "outputs")
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures")

# Human-readable labels, ordered the way the panel is reported.
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


def ensure_dirs() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FIGURE_DIR, exist_ok=True)


def connect(db_path: Optional[str] = None) -> sqlite3.Connection:
    path = db_path or DB_PATH
    if not os.path.exists(path):
        raise FileNotFoundError(
            "Database not found at %s. Run `python load_data.py` (or `make pipeline`) first." % path
        )
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def query(sql: str, params=(), db_path: Optional[str] = None) -> pd.DataFrame:
    conn = connect(db_path)
    try:
        return pd.read_sql_query(sql, conn, params=params)
    finally:
        conn.close()


def analysis_frame(db_path: Optional[str] = None) -> pd.DataFrame:
    """The full denormalised sample x population table (one row per measurement)."""
    df = query("SELECT * FROM v_analysis_base", db_path=db_path)
    df["population_label"] = df["population"].map(lambda p: POPULATION_LABELS.get(p, p))
    return df


def population_order(df: pd.DataFrame) -> list:
    """Populations in panel order, restricted to what is actually present."""
    present = list(dict.fromkeys(df["population"]))
    ordered = [p for p in POPULATION_LABELS if p in present]
    return ordered + [p for p in present if p not in ordered]


def write_csv(df: pd.DataFrame, name: str) -> str:
    ensure_dirs()
    path = os.path.join(OUTPUT_DIR, name)
    df.to_csv(path, index=False)
    return path


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

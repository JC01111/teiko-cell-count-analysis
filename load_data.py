#!/usr/bin/env python3
"""Initialise the SQLite database and load cell-count.csv into it.

Usage
-----
    python load_data.py                       # data/cell-count.csv -> db/cell_counts.db
    python load_data.py --csv other.csv --db /tmp/foo.db
    python load_data.py --append              # keep existing rows, upsert new ones

The loader is *column-driven*: any CSV column that is not a known clinical
attribute is treated as a cell population and registered in `cell_population`.
Dropping a new marker column into the CSV therefore needs no code or schema
change.
"""
from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import sys
from typing import Dict, List, Optional

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(ROOT, "data", "cell-count.csv")
DEFAULT_DB = os.path.join(ROOT, "db", "cell_counts.db")
SCHEMA_SQL = os.path.join(ROOT, "db", "schema.sql")

# Columns that describe the patient/sample rather than a cell population.
METADATA_COLUMNS = {
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
}

DISPLAY_NAMES = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


def _clean(value: Optional[str]) -> Optional[str]:
    """Trim whitespace; treat empty strings / NA markers as SQL NULL."""
    if value is None:
        return None
    value = value.strip()
    # NB: "none" is NOT treated as missing -- it is a real treatment value
    # ("no treatment") in this dataset.
    if value == "" or value.lower() in {"na", "n/a", "nan", "null"}:
        return None
    return value


def _to_int(value: Optional[str], field: str, line: int) -> Optional[int]:
    value = _clean(value)
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        raise ValueError("line %d: %r is not a valid integer for %s" % (line, value, field))


class _Cache:
    """Insert-or-fetch helper for the small dimension tables."""

    def __init__(self, conn: sqlite3.Connection, table: str, key_col: str, id_col: str):
        self.conn, self.table, self.key_col, self.id_col = conn, table, key_col, id_col
        self.map: Dict[str, int] = {}
        for row in conn.execute("SELECT %s, %s FROM %s" % (id_col, key_col, table)):
            self.map[row[1]] = row[0]

    def get(self, key: Optional[str], **extra) -> Optional[int]:
        if key is None:
            return None
        if key in self.map:
            return self.map[key]
        cols = [self.key_col] + list(extra)
        sql = "INSERT INTO %s (%s) VALUES (%s)" % (
            self.table,
            ", ".join(cols),
            ", ".join("?" for _ in cols),
        )
        cur = self.conn.execute(sql, [key] + list(extra.values()))
        self.map[key] = cur.lastrowid
        return cur.lastrowid


def create_schema(conn: sqlite3.Connection) -> None:
    with open(SCHEMA_SQL, "r", encoding="utf-8") as fh:
        conn.executescript(fh.read())


def load_csv(conn: sqlite3.Connection, csv_path: str) -> Dict[str, int]:
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("%s appears to be empty" % csv_path)

        header = [h.strip() for h in reader.fieldnames]
        missing = METADATA_COLUMNS - set(header)
        if missing:
            raise ValueError("missing required column(s): %s" % ", ".join(sorted(missing)))
        population_cols: List[str] = [h for h in header if h not in METADATA_COLUMNS]
        if not population_cols:
            raise ValueError("no cell-population columns found in %s" % csv_path)

        conditions = _Cache(conn, "condition", "condition_name", "condition_id")
        treatments = _Cache(conn, "treatment", "treatment_name", "treatment_id")
        types = _Cache(conn, "sample_type", "sample_type_name", "sample_type_id")
        projects = _Cache(conn, "project", "project_code", "project_id")
        subjects = _Cache(conn, "subject", "subject_code", "subject_id")
        populations = _Cache(conn, "cell_population", "population_name", "population_id")

        pop_ids = {}
        for order, col in enumerate(population_cols):
            pop_ids[col] = populations.get(
                col, display_name=DISPLAY_NAMES.get(col, col.replace("_", " ")), sort_order=order
            )

        enrollments: Dict[tuple, int] = {}
        for row in conn.execute(
            "SELECT enrollment_id, subject_id, project_id, condition_id, treatment_id FROM enrollment"
        ):
            enrollments[row[1:]] = row[0]

        n_samples = n_counts = 0
        for line, raw in enumerate(reader, start=2):
            rec = {k.strip(): v for k, v in raw.items() if k is not None}
            subject_code = _clean(rec["subject"])
            sample_code = _clean(rec["sample"])
            if subject_code is None or sample_code is None:
                raise ValueError("line %d: subject and sample are required" % line)

            sex = _clean(rec["sex"])
            subject_id = subjects.get(subject_code, sex=sex)
            project_id = projects.get(_clean(rec["project"]))
            condition_id = conditions.get(_clean(rec["condition"]))
            treatment_id = treatments.get(_clean(rec["treatment"]))
            key = (subject_id, project_id, condition_id, treatment_id)
            enrollment_id = enrollments.get(key)
            if enrollment_id is None:
                cur = conn.execute(
                    "INSERT INTO enrollment (subject_id, project_id, condition_id,"
                    " treatment_id, response, age_at_enrollment) VALUES (?,?,?,?,?,?)",
                    key + (_clean(rec["response"]), _to_int(rec["age"], "age", line)),
                )
                enrollment_id = cur.lastrowid
                enrollments[key] = enrollment_id

            cur = conn.execute(
                "INSERT INTO sample (sample_code, enrollment_id, sample_type_id,"
                " time_from_treatment_start) VALUES (?,?,?,?)"
                " ON CONFLICT(sample_code) DO UPDATE SET"
                " enrollment_id=excluded.enrollment_id,"
                " sample_type_id=excluded.sample_type_id,"
                " time_from_treatment_start=excluded.time_from_treatment_start",
                (
                    sample_code,
                    enrollment_id,
                    types.get(_clean(rec["sample_type"])),
                    _to_int(rec["time_from_treatment_start"], "time_from_treatment_start", line),
                ),
            )
            sample_id = cur.lastrowid or conn.execute(
                "SELECT sample_id FROM sample WHERE sample_code = ?", (sample_code,)
            ).fetchone()[0]
            n_samples += 1

            payload = []
            for col in population_cols:
                count = _to_int(rec.get(col), col, line)
                if count is None:
                    continue  # marker not measured on this panel -> simply no row
                payload.append((sample_id, pop_ids[col], count))
            conn.executemany(
                "INSERT INTO sample_cell_count (sample_id, population_id, cell_count)"
                " VALUES (?,?,?) ON CONFLICT(sample_id, population_id)"
                " DO UPDATE SET cell_count = excluded.cell_count",
                payload,
            )
            n_counts += len(payload)

    return {"samples": n_samples, "cell_counts": n_counts, "populations": len(population_cols)}


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", default=DEFAULT_CSV, help="path to cell-count.csv")
    parser.add_argument("--db", default=DEFAULT_DB, help="path to the SQLite database file")
    parser.add_argument(
        "--append",
        action="store_true",
        help="load into an existing database instead of rebuilding it from scratch",
    )
    args = parser.parse_args(argv)

    csv_path = args.csv
    if not os.path.exists(csv_path):
        fallback = os.path.join(ROOT, "cell-count.csv")
        if os.path.exists(fallback):
            csv_path = fallback
        else:
            parser.error("CSV not found: %s" % args.csv)

    os.makedirs(os.path.dirname(os.path.abspath(args.db)) or ".", exist_ok=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        with conn:
            if not args.append:
                create_schema(conn)
            stats = load_csv(conn, csv_path)
        conn.execute("ANALYZE")
    finally:
        conn.close()

    print("Loaded %s -> %s" % (csv_path, args.db))
    print(
        "  %(samples)d samples, %(populations)d cell populations, %(cell_counts)d count rows"
        % stats
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

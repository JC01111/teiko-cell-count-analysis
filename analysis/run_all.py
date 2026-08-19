#!/usr/bin/env python3
"""Run the whole analysis pipeline (Parts 2-4) and write every output."""
from __future__ import annotations

import sys

from analysis import frequencies, statistics_report, subsets


def main() -> int:
    for step in (frequencies.main, statistics_report.main, subsets.main):
        sys.argv = [sys.argv[0]]
        rc = step()
        if rc:
            return rc
    print("\nPipeline complete. See outputs/ for tables and outputs/figures/ for plots.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

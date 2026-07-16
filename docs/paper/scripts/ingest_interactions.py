#!/usr/bin/env python3
"""Extract cuda_interactions_test residuals from a SLURM cuda_test_*.out file.

Captures the 5 sub-tests: NC+CC, tau regen, Glashow, E^-2, E^-3.
Output: CSV mirroring data/interactions_*.csv schema.

Usage:
    python ingest_interactions.py path/to/cuda_test_<jobid>.out [--commit ...]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

TEST_HEADER = re.compile(r"^--- Test (\d+): (.+?) ---")
ABS_LINE = re.compile(r"Max absolute error: ([\d.eE+-]+)")
REL_LINE = re.compile(r"Max relative error: ([\d.eE+-]+)")
RESULT = re.compile(r"Result: (PASS|FAIL)")


def parse(path: Path) -> list[dict]:
    rows: list[dict] = []
    in_interactions = False
    cur: dict | None = None
    for line in path.read_text().splitlines():
        if "=== Running cuda_interactions_test ===" in line:
            in_interactions = True
            continue
        if not in_interactions:
            continue
        if "=== Running benchmark_cuda ===" in line:
            break
        m = TEST_HEADER.match(line)
        if m and "GPU vs CPU" in m.group(2):
            if cur:
                rows.append(cur)
            cur = {"test": int(m.group(1)), "description": m.group(2).strip()}
            continue
        if cur is None:
            continue
        if (m := ABS_LINE.search(line)):
            cur["max_abs_err"] = float(m.group(1))
        elif (m := REL_LINE.search(line)):
            cur["max_rel_err"] = float(m.group(1))
        elif (m := RESULT.search(line)):
            cur["status"] = m.group(1)
    if cur and "status" in cur:
        rows.append(cur)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", type=Path)
    ap.add_argument("--commit", default="UNKNOWN")
    args = ap.parse_args()

    rows = parse(args.infile)
    if not rows:
        print("no interaction-test rows found", file=sys.stderr)
        return 2

    cols = ["test", "description", "max_abs_err", "max_rel_err", "status", "commit"]
    print(",".join(cols))
    for r in rows:
        r["commit"] = args.commit
        print(",".join(str(r.get(c, "")) for c in cols))
    return 0


if __name__ == "__main__":
    sys.exit(main())

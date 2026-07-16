#!/usr/bin/env python3
"""Extract the benchmark_cuda sweep table from a SLURM cuda_test_*.out file.

Usage:
    python ingest_benchmark.py path/to/cuda_test_<jobid>.out [--gpu-name "..."] [--partition "..."] > out.csv

The benchmark table looks like (after `=== Running benchmark_cuda ===`):

   ncz    ne      CPU (ms)      GPU (ms)   Speedup   Max Rel Err  Status
------------------------------------------------------------------------
     5    20           8.8         166.7      0.05x      6.91e-03    PASS
    ...

Output: CSV columns matching data/benchmark_*.csv schema.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROW = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)x\s+([\d.eE+-]+)\s+(PASS|FAIL)"
)


def parse(path: Path) -> list[dict]:
    rows = []
    in_bench = False
    for line in path.read_text().splitlines():
        if "=== Running benchmark_cuda ===" in line:
            in_bench = True
            continue
        if not in_bench:
            continue
        m = ROW.match(line)
        if not m:
            continue
        ncz, ne, cpu, gpu, sp, err, status = m.groups()
        rows.append(
            {
                "ncz": int(ncz),
                "ne": int(ne),
                "points": int(ncz) * int(ne),
                "cpu_ms": float(cpu),
                "gpu_ms": float(gpu),
                "speedup": float(sp),
                "max_rel_err": float(err),
                "status": status,
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", type=Path)
    ap.add_argument("--gpu-name", default="UNKNOWN")
    ap.add_argument("--partition", default="UNKNOWN")
    ap.add_argument("--commit", default="UNKNOWN")
    args = ap.parse_args()

    if not args.infile.exists():
        print(f"no such file: {args.infile}", file=sys.stderr)
        return 1

    jobid_match = re.search(r"cuda_test(?:_full_a100)?_(\d+)\.out", args.infile.name)
    jobid = jobid_match.group(1) if jobid_match else "UNKNOWN"

    rows = parse(args.infile)
    if not rows:
        print(f"no benchmark rows found in {args.infile}", file=sys.stderr)
        return 2

    cols = [
        "gpu",
        "partition",
        "ncz",
        "ne",
        "points",
        "cpu_ms",
        "gpu_ms",
        "speedup",
        "max_rel_err",
        "status",
        "slurm_jobid",
        "commit",
    ]
    print(",".join(cols))
    for r in rows:
        r["gpu"] = args.gpu_name
        r["partition"] = args.partition
        r["slurm_jobid"] = jobid
        r["commit"] = args.commit
        print(",".join(str(r[c]) for c in cols))

    return 0


if __name__ == "__main__":
    sys.exit(main())

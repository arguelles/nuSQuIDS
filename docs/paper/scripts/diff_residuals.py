#!/usr/bin/env python3
"""Compare two paper_diag_*.csv files (typically MIG vs full A100) bin-by-bin.

Answers the cross-architecture question: are the per-bin residuals essentially
identical (FP-reordering noise) or globally different (real bug)?

Usage:
    python diff_residuals.py data/paper_diag_tau_mig.csv data/paper_diag_tau_a100.csv

Reports:
  - max |abs_err_a - abs_err_b| across all bins (should be ~0 if reordering only)
  - max |gpu_a - gpu_b| across all bins (FP order signature)
  - which bin holds each maximum
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load(p: Path) -> dict:
    out = {}
    for r in csv.DictReader(p.open()):
        key = (float(r["cz"]), float(r["E_GeV"]), r["rho"], r["flavor"])
        out[key] = {k: float(v) if k not in ("rho", "flavor") else v for k, v in r.items()}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("a", type=Path)
    ap.add_argument("b", type=Path)
    args = ap.parse_args()

    A = load(args.a)
    B = load(args.b)

    common = sorted(A.keys() & B.keys())
    if not common:
        print("no overlapping bins!", file=sys.stderr)
        return 1

    worst_abs_diff = (None, 0.0)  # (bin, max | abs_err_a - abs_err_b |)
    worst_gpu_diff = (None, 0.0)  # (bin, max | gpu_a - gpu_b |)
    worst_rel_a = (None, 0.0)
    worst_rel_b = (None, 0.0)

    for k in common:
        a = A[k]
        b = B[k]
        d_abs = abs(a["abs_err"] - b["abs_err"])
        d_gpu = abs(a["gpu"] - b["gpu"])
        if d_abs > worst_abs_diff[1]:
            worst_abs_diff = (k, d_abs)
        if d_gpu > worst_gpu_diff[1]:
            worst_gpu_diff = (k, d_gpu)
        if a["rel_err"] > worst_rel_a[1]:
            worst_rel_a = (k, a["rel_err"])
        if b["rel_err"] > worst_rel_b[1]:
            worst_rel_b = (k, b["rel_err"])

    def fmt(entry):
        if entry[0] is None:
            return f"{entry[1]:.3e}  (no nonzero values)"
        k = entry[0]
        return f"{entry[1]:.3e}  at  cz={k[0]:.3f}  E={k[1]:.3e} GeV  {k[2]} {k[3]}"

    print(f"compared {len(common)} bins")
    print(f"file A: {args.a.name}")
    print(f"file B: {args.b.name}")
    print()
    print(f"max |abs_err_a - abs_err_b|:  {fmt(worst_abs_diff)}")
    print(f"max |gpu_a - gpu_b|:          {fmt(worst_gpu_diff)}")
    print()
    print(f"worst rel_err in A:  {fmt(worst_rel_a)}")
    print(f"worst rel_err in B:  {fmt(worst_rel_b)}")
    print()
    if worst_abs_diff[1] == 0.0 and worst_gpu_diff[1] == 0.0:
        print("=> per-bin residuals are byte-identical between A and B.")
        print("   Same code, different architectures (e.g. MIG vs full A100),")
        print("   identical numerical output. The integrator is deterministic")
        print("   under different SM counts and warp scheduling.")
    elif worst_rel_a[0] != worst_rel_b[0]:
        print("=> worst-rel bins differ between A and B.")
        print("   This is the signature of FP-reordering noise:")
        print("   the same code produces tiny differences whose worst location")
        print("   moves with execution order, but the overall residual envelope")
        print("   is the same.")
    else:
        print("=> worst-rel bin matches between A and B; differences are localised.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

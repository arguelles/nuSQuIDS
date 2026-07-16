#!/usr/bin/env python3
"""Render the cuda_interactions_test residuals as a small TeX table or matplotlib bar chart.

Usage:
    python plot_residuals.py data/interactions_3bab4cc.csv -o figures/residuals.pdf
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("figures/residuals.pdf"))
    args = ap.parse_args()

    rows = list(csv.DictReader(args.infile.open()))
    labels = [f"{r['test']}: {r['description']}" for r in rows]
    abs_err = [float(r["max_abs_err"]) for r in rows]
    rel_err = [float(r["max_rel_err"]) for r in rows]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = range(len(labels))
    ax.bar([xi - 0.2 for xi in x], abs_err, width=0.4, label="max abs err")
    ax.bar([xi + 0.2 for xi in x], rel_err, width=0.4, label="max rel err")
    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.axhline(5e-2, linestyle="--", color="C2", label="rel tol (Tests 2-5)")
    ax.axhline(1e-2, linestyle="--", color="C3", label="abs tol (Tests 2-5)")
    ax.set_ylabel("residual")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

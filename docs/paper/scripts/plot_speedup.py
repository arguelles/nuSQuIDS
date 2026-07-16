#!/usr/bin/env python3
"""Plot GPU/CPU speedup vs grid size from one or more benchmark CSVs.

Usage:
    python plot_speedup.py data/benchmark_12038685.csv [data/benchmark_otherjob.csv ...] -o figures/speedup.pdf

X axis: ncz * ne (total grid points), log scale.
Y axis: speedup (CPU ms / GPU ms), log scale.
Each input CSV becomes a series labeled by its `gpu` column.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load(path: Path) -> tuple[str, list[int], list[float]]:
    points = []
    speedups = []
    label = path.stem
    with path.open() as f:
        for row in csv.DictReader(f):
            if row["status"] != "PASS":
                continue
            label = row.get("gpu") or label
            points.append(int(row["points"]))
            speedups.append(float(row["speedup"]))
    return label, points, speedups


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("figures/speedup.pdf"))
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(7, 5))
    for p in args.inputs:
        label, pts, sp = load(p)
        ax.loglog(pts, sp, marker="o", linestyle="-", label=label)

    ax.axhline(1.0, color="k", linestyle="--", linewidth=0.8, label="GPU = CPU")
    ax.set_xlabel(r"Grid points $n_{c_\theta} \times n_E$")
    ax.set_ylabel(r"Speedup (CPU time / GPU time)")
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

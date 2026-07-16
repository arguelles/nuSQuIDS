#!/usr/bin/env python3
"""Plot CPU vs GPU per-bin agreement from paper_diag_*.csv files.

Two outputs:
  1. heatmap of relative error in (cosθ, log10 E) at fixed flavor/rho
  2. scatter of GPU vs CPU values across all bins, with diagonal

Usage:
    python plot_perbin.py data/paper_diag_tau_<jobid>.csv -o figures/perbin_tau.pdf

The CSV is expected to have columns:
  cz, E_GeV, rho, flavor, cpu, gpu, abs_err, rel_err
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load(path: Path):
    rows = list(csv.DictReader(path.open()))
    czs = sorted(set(float(r["cz"]) for r in rows))
    es = sorted(set(float(r["E_GeV"]) for r in rows))
    flavors = sorted(set(r["flavor"] for r in rows))
    rhos = sorted(set(r["rho"] for r in rows))
    return rows, czs, es, flavors, rhos


def heatmap(ax, rows, flv, rho, czs, es, value_key):
    data = np.full((len(czs), len(es)), np.nan)
    for r in rows:
        if r["flavor"] != flv or r["rho"] != rho:
            continue
        i = czs.index(float(r["cz"]))
        j = es.index(float(r["E_GeV"]))
        data[i, j] = float(r[value_key])
    if value_key == "rel_err":
        data = np.where(data > 1e-15, data, 1e-15)
        im = ax.pcolormesh(np.log10(es), czs, data,
                           shading="auto", norm=__import__("matplotlib").colors.LogNorm(vmin=1e-6, vmax=1.0),
                           cmap="viridis")
    else:
        im = ax.pcolormesh(np.log10(es), czs, data, shading="auto", cmap="viridis")
    ax.set_xlabel(r"$\log_{10}(E / \mathrm{GeV})$")
    ax.set_ylabel(r"$\cos\theta$")
    ax.set_title(f"{flv} ({rho})")
    return im


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("figures/perbin.pdf"))
    ap.add_argument("--flavor", default="nu_mu", help="flavor for heatmap (default nu_mu)")
    args = ap.parse_args()

    rows, czs, es, flavors, rhos = load(args.infile)

    # 2x1 grid: heatmap of rel_err for the chosen flavor (nu and nubar)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    last_im = None
    for ax, rho in zip(axes, ["nu", "nubar"]):
        last_im = heatmap(ax, rows, args.flavor, rho, czs, es, "rel_err")
    cbar = fig.colorbar(last_im, ax=axes, label="|GPU - CPU| / |CPU|")
    cbar.ax.set_ylabel("relative error", rotation=90)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out)
    print(f"wrote {args.out}")

    # Scatter plot, separate file
    scat_path = args.out.with_name(args.out.stem + "_scatter.pdf")
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    cpu = np.array([float(r["cpu"]) for r in rows])
    gpu = np.array([float(r["gpu"]) for r in rows])
    mask = (np.abs(cpu) > 1e-12) & (np.abs(gpu) > 1e-12)
    ax.loglog(np.abs(cpu[mask]), np.abs(gpu[mask]), '.', alpha=0.3, markersize=2)
    lo = max(np.abs(cpu[mask]).min(), 1e-12)
    hi = max(np.abs(cpu[mask]).max(), 1.0)
    ax.loglog([lo, hi], [lo, hi], 'r-', linewidth=0.8, label="GPU = CPU")
    ax.set_xlabel(r"|CPU EvalFlavor|")
    ax.set_ylabel(r"|GPU EvalFlavor|")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(scat_path)
    print(f"wrote {scat_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

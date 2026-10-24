#!/usr/bin/env python3
"""Work-precision diagrams from a sweep CSV (one panel per workload).

usage: plot_wp.py <sweep.csv> <out.pdf>
"""
import csv
import math
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (method, label, color, linestyle, marker)
SHOW = [
    ("rkf45", "GSL rkf45 (default)", "black", "-", "o"),
    ("rkck", "GSL rkck", "0.55", "--", "s"),
    ("rk8pd", "GSL rk8pd", "0.55", ":", "^"),
    ("n_rkf45:I", "native RKF45 (same method, new driver)", "#1f77b4", "--", "o"),
    ("n_dp5:I", "native DP5", "#2ca02c", "-", "s"),
    ("n_vern7:I", "native Vern7", "#d62728", "-", "D"),
    ("n_vern9:I", "native Vern9", "#9467bd", "-", "v"),
    ("ni_vern7:I", "per-node Vern7", "#d62728", ":", "D"),
    ("ni_magnus4", "per-node Magnus-4", "#8c564b", "-.", "x"),
    ("ni_modmagnus", "per-node modified Magnus", "#ff7f0e", "-", "*"),
]
TITLES = {
    "earth_osc": "core-crossing, 0.5-500 GeV",
    "const_osc": "constant density, 0.5-500 GeV",
    "atm_osc": "atmospheric grid (20 zeniths), 1-1000 GeV",
    "earth_int": "interactions + tau regen., 0.1-10^4 TeV",
    "sterile": "3+1 sterile, 0.1-100 TeV",
    "earth_subgev": "core-crossing, 0.1-2 GeV",
    "earth_mev": "core-crossing, 2-20 MeV",
}

rows = list(csv.DictReader(open(sys.argv[1])))
pts = defaultdict(lambda: defaultdict(list))
for r in rows:
    try:
        e, t = float(r["maxerr"]), float(r["time_s"])
    except ValueError:
        continue
    if not (math.isnan(e) or math.isnan(t)) and e > 0:
        pts[r["workload"]][r["method"]].append((e, t))

wls = [w for w in TITLES if w in pts]
ncol = 2
nrow = (len(wls) + ncol - 1)//ncol
fig, axes = plt.subplots(nrow, ncol, figsize=(10, 3.6*nrow), squeeze=False)
for ax, w in zip(axes.flat, wls):
    for m, label, color, ls, mk in SHOW:
        p = sorted(pts[w].get(m, []))
        if not p:
            continue
        ax.plot([e for e, _ in p], [t for _, t in p], color=color, ls=ls, marker=mk,
                ms=4, lw=1.2, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(TITLES[w], fontsize=10)
    ax.set_xlabel(r"max $|\Delta P|$ against converged reference")
    ax.set_ylabel("time in EvolveState [s]")
    ax.grid(True, which="both", alpha=0.25)
    ax.invert_xaxis()
for ax in list(axes.flat)[len(wls):]:
    ax.axis("off")
handles, labels = axes.flat[0].get_legend_handles_labels()
for ax in axes.flat:
    h, l = ax.get_legend_handles_labels()
    for hh, ll in zip(h, l):
        if ll not in labels:
            handles.append(hh)
            labels.append(ll)
fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(sys.argv[2])
print("wrote", sys.argv[2])

#!/usr/bin/env python3
"""Work-precision frontier from sweep CSVs.

For each workload and target accuracy, report every method's cheapest run that
achieved max error <= target (at whatever tolerance that took), and its
speedup over the baseline method.

usage: frontier.py <sweep.csv> [more.csv ...] [--baseline rkf45] [--targets 1e-4,1e-6,1e-8]
       [--metric time_s|nfev]
"""
import argparse
import csv
import math
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("csv", nargs="+")
ap.add_argument("--baseline", default="rkf45")
ap.add_argument("--targets", default="1e-4,1e-5,1e-6,1e-7,1e-8,1e-9")
ap.add_argument("--metric", default="time_s")
ap.add_argument("--top", type=int, default=6)
args = ap.parse_args()

rows = []
for f in args.csv:
    rows += list(csv.DictReader(open(f)))
targets = [float(x) for x in args.targets.split(",")]

data = defaultdict(lambda: defaultdict(list))  # wl -> method -> [(err, cost)]
for r in rows:
    try:
        err, cost = float(r["maxerr"]), float(r[args.metric])
    except ValueError:
        continue
    if math.isnan(err) or math.isnan(cost):
        continue
    data[r["workload"]][r["method"]].append((err, cost))


def best(points, target):
    ok = [c for e, c in points if e <= target]
    return min(ok) if ok else None


for wl in data:
    print("== %s (metric: %s, baseline: %s)" % (wl, args.metric, args.baseline))
    print("%-9s %-12s %10s   %s" % ("target", "baseline", "", "fastest methods (speedup vs baseline)"))
    for t in targets:
        base = best(data[wl].get(args.baseline, []), t)
        ranking = []
        for m, pts in data[wl].items():
            c = best(pts, t)
            if c is not None:
                ranking.append((c, m))
        ranking.sort()
        cells = []
        for c, m in ranking[: args.top]:
            sp = ("%.2fx" % (base / c)) if base else "  -  "
            cells.append("%s %s" % (m, sp))
        print("%-9.0e %-12s %10s   %s" % (t, "%.4g" % base if base else "n/a", "", ", ".join(cells)))
    print()

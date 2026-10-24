#!/usr/bin/env python3
"""Rank methods by geometric-mean speedup at equal achieved accuracy.

For every (workload, target accuracy) the cost of a method is its cheapest run
with max error <= target. Speedups are taken against (a) the default GSL rkf45
and (b) the best GSL stepper at that point. Only points every compared method
reached are used, so the averages are over the same set.

usage: rank.py <sweep.csv> [--targets 1e-5,1e-6,1e-7,1e-8] [--workloads a,b]
"""
import argparse
import csv
import math
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("csv", nargs="+")
ap.add_argument("--targets", default="1e-5,1e-6,1e-7,1e-8")
ap.add_argument("--workloads", default="")
ap.add_argument("--metric", default="time_s")
args = ap.parse_args()

GSL = {"rkf45", "rkck", "rk8pd", "rk4", "msadams"}
targets = [float(x) for x in args.targets.split(",")]
data = defaultdict(lambda: defaultdict(list))
for f in args.csv:
    for r in csv.DictReader(open(f)):
        try:
            e, c = float(r["maxerr"]), float(r[args.metric])
        except ValueError:
            continue
        if not (math.isnan(e) or math.isnan(c)):
            data[r["workload"]][r["method"]].append((e, c))
wls = [w for w in data if not args.workloads or w in args.workloads.split(",")]


def best(pts, t):
    ok = [c for e, c in pts if e <= t]
    return min(ok) if ok else None


methods = sorted({m for w in wls for m in data[w]})
score = {}
for m in methods:
    sp_def, sp_gsl, n = [], [], 0
    for w in wls:
        for t in targets:
            c = best(data[w].get(m, []), t)
            base = best(data[w].get("rkf45", []), t)
            gsl = [best(data[w][g], t) for g in GSL if g in data[w]]
            gsl = [x for x in gsl if x is not None]
            if c is None or base is None or not gsl:
                continue
            sp_def.append(math.log(base / c))
            sp_gsl.append(math.log(min(gsl) / c))
            n += 1
    if n:
        score[m] = (math.exp(sum(sp_def) / n), math.exp(sum(sp_gsl) / n), n)

total = len(wls) * len(targets)
print("workloads: %s; targets: %s" % (",".join(wls), args.targets))
print("%-16s %10s %12s %8s" % ("method", "vs rkf45", "vs best GSL", "points"))
for m, (a, b, n) in sorted(score.items(), key=lambda kv: -kv[1][0]):
    flag = "" if n == total else "  (reached %d/%d)" % (n, total)
    print("%-16s %9.2fx %11.2fx %8d%s" % (m, a, b, n, flag))

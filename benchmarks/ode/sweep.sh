#!/bin/bash
# Work-precision sweep: every (workload, method, tolerance) combination, with the
# max absolute error against the converged reference in ref/.
#
# usage: sweep.sh <out.csv> <reps> "<workloads>" "<methods>" "<tolerances>"
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="$1"; reps="$2"; workloads="$3"; methods="$4"; tols="$5"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
[ -f "$out" ] || echo "workload,method,tol,time_s,nfev,steps,rejected,maxerr" > "$out"
for w in $workloads; do
  ref="$here/ref/${w}_rk8pd_1e-14.txt"
  for m in $methods; do
    for t in $tols; do
      line="$("${BENCH:-$here/ode_bench}" "$w" "$m" "$t" "$tmp/r.txt" "$reps" 2>&1 | tail -1)"
      if [[ "$line" == *FAILED* || "$line" != *time_s=* ]]; then
        echo "$w,$m,$t,NaN,NaN,NaN,NaN,NaN" >> "$out"
        continue
      fi
      time_s="$(sed -E 's/.*time_s=([^ ]+).*/\1/' <<< "$line")"
      nfev="$(sed -E 's/.*nfev=([^ ]+).*/\1/' <<< "$line")"
      steps="$(sed -E 's/.*steps=([^ ]+).*/\1/' <<< "$line")"
      rej="$(sed -E 's/.*rejected=([^ ]+).*/\1/' <<< "$line")"
      err="$(python3 "$here/maxerr.py" "$tmp/r.txt" "$ref")"
      echo "$w,$m,$t,$time_s,$nfev,$steps,$rej,$err" >> "$out"
    done
  done
done

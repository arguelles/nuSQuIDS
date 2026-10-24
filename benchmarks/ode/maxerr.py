#!/usr/bin/env python3
"""Max absolute difference between two ode_bench output files."""
import sys

a = [float(x) for x in open(sys.argv[1])]
b = [float(x) for x in open(sys.argv[2])]
if len(a) != len(b):
    sys.exit("length mismatch: %d vs %d" % (len(a), len(b)))
print("%.3e" % max(abs(x - y) for x, y in zip(a, b)))

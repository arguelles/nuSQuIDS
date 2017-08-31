#!/usr/bin/python

import numpy as np
import sys


if len(sys.argv)<2:
    raise Exception("An imput file needed")
fl=sys.argv[1]

file=open(fl, 'r') 
for f in file:
    if f[0]=='a':
        cz=0.5*(float(f[23:28])+float(f[32:37]))
    elif f[:2]!=' E':
        print cz,f[:-1]

file.close()

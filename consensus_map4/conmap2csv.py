#!/usr/bin/python

'''
combine individual consensus LGs into a single csv file map
'''

import sys
import os

outfile = sys.argv[1]
inpfiles = sys.argv[2:]

data = []

for inp in inpfiles:
    base = os.path.basename(inp)
    lg = os.path.splitext(base)[0]
    f = open(inp)
    f.readline()
    for line in f:
        tok = line.strip().split(',')
        uid = tok[0]
        pos = float(tok[1])
        data.append([uid,lg,pos])
    f.close()

data.sort(key=lambda x:x[2])
data.sort(key=lambda x:x[1])

fout = open(outfile,'wb')
for row in data: fout.write('%s,%s,%.4f\n'%(row[0],row[1],row[2]))
fout.close()

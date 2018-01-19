#!/usr/bin/python

'''
count missed calls per row (sample or marker depending on the input file)
'''

import sys

inp   = sys.argv[1]
batch = sys.argv[2]

f = open(inp)
header = f.readline() #header
ncols = len(header.strip().split('\t')) - 1

for line in f:
    tok = line.strip().split('\t')
    assert len(tok) == ncols + 1
    name = tok[0]
    missing = [int(x) for x in tok[1:]].count(-1)
    print name,batch,missing,ncols
f.close()

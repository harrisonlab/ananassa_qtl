#!/usr/bin/python

#
# get just the two parental genotypes
# convert from probe to snp id
#

import sys

inp = sys.argv[1]
p2sfile="/home/vicker/istraw90/IStraw90.r1.ps2snp_map.ps.fixed"

p2s = {}
f = open(p2sfile)
f.readline()
for line in f:
    tok = line.strip().split()
    pid = tok[0][3:]
    sid = tok[1]
    p2s[pid] = sid
f.close()

f = open(inp)
f.readline()
for line in f:
    tok = line.strip().split('\t')
    pid = tok[0].split('-')[1]
    sid = p2s[pid]
    mat = tok[1]
    pat = tok[2]
    print '\t'.join([sid,mat,pat])
f.close()

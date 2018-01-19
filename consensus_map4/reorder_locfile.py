#!/usr/bin/python

#
# reorder loc file according to a given map file
# ignore any markers in the map file not in the loc file
# drop markers in the loc file not in the map file
# locfile should correspond to a single linkage group only
#

import sys

inpfile=sys.argv[1]
orderfile=sys.argv[2]
target_lg=sys.argv[3]
outfile=sys.argv[4]

mid2pos = {}

f = open(orderfile)
for line in f:
    tok = line.strip().split(',')
    mid = tok[0]
    lg = tok[1]
    pos = float(tok[2])

    if lg != target_lg or target_lg == "-": continue #ignore this marker from wrong lg

    mid2pos[mid] = pos

f.close()

data = []
f = open(inpfile)
for line in f:
    mid = line.split()[0]
    if not mid in mid2pos: continue #drop markers without a position in the order file
    data.append([mid,line])
f.close()

#sort by marker position
data.sort(key=lambda x:mid2pos[x[0]])

fout = open(outfile,'wb')
for row in data: fout.write(row[1])
fout.close()

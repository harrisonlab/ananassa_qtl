#!/usr/bin/python

#
# compare rf / LOD values of markers in a linkage group with grouping in a reference map
# report markers which seem to be spanning two genuine lgs
# report only markers spanning homeologous lgs
#

import sys

rflodfile = sys.argv[1] #rflod for a single linkage groups: 
mapfile = sys.argv[2]   #map data from an existing map

#load linkage group designation of markers in existing map
mlg = {}
f = open(mapfile)
for line in f:
    tok = line.strip().split(',')
    mlg[tok[0]] = tok[1]
f.close()

def update_item(a,b,mlg,mscore,lod):
    if not a in mlg: return
    
    lg = mlg[a]
    
    if not lg in mscore[b]: mscore[b][lg] = [0.0, 0]
    
    mscore[b][lg][0] += lod
    mscore[b][lg][1] += 1

#lod rflod data, record max lod of each marker to each LG
mscore = {}
f = open(rflodfile)
for line in f:
    tok = line.strip().split()
    m1 = tok[0]
    m2 = tok[1]
    rf = float(tok[2])
    lod = float(tok[3])

    if not m1 in mscore: mscore[m1] = {}
    if not m2 in mscore: mscore[m2] = {}
    
    update_item(m1,m2,mlg,mscore,lod)
    update_item(m2,m1,mlg,mscore,lod)
    
f.close()

for x in mscore:
    if len(mscore[x]) == 1: continue #ignore if only linked to one linkage group
    keys = mscore[x].keys()
    for i,lg1 in enumerate(keys):
        for lg2 in keys[i+1:]:
            #only report linkage to pairs of homeologous linkage groups
            #if lg1[:-1] != lg2[:-1]: continue
            
            print x,lg1,mscore[x][lg1][0],mscore[x][lg1][1],lg2,mscore[x][lg2][0],mscore[x][lg2][1]

#!/usr/bin/python

#
# merge split groups which show sufficient linkage
#

#../tmpmergegrps/${lg}.rflod ${outdir} *.loc

import sys
import os
from unionfind import UnionFind

rflod = sys.argv[1]
inpdir = sys.argv[2]
outdir = sys.argv[3]
lg = sys.argv[4]
frags = sys.argv[5].split()

#load lod values, retain only the highest frag-to-frag lod
ff2lod = {}
f = open(rflod)
for line in f:
    tok = line.strip().split()
    frag1 = tok[0].split('-')[0]
    frag2 = tok[1].split('-')[0]
    if frag1 > frag2: frag1,frag2=frag2,frag1
    key = frag1+':'+frag2
    lod = float(tok[3])

    if not key in ff2lod:   ff2lod[key] = lod
    elif lod > ff2lod[key]: ff2lod[key] = lod

f.close()

#union-find frags into groups
uf = UnionFind()
for key in ff2lod:
    frag1,frag2 = key.split(':')
    uf.union(frag1,frag2)

merged = {}
for frag in frags:
    mid = uf[frag]
    if not mid in merged: merged[mid] = []
    merged[mid].append(frag)

for mid in merged:
    outfile = outdir + '/' + lg + '.' + '_'.join(merged[mid]) + '.loc'
    inplist = ['%s/%s.loc'%(inpdir,x) for x in merged[mid]]

    cmd = 'cat ' + ' '.join(inplist) + ' > ' + outfile
    assert os.system(cmd) == 0


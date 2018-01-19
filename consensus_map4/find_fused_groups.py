#!/usr/bin/python

#
# work out which provisional LGs are a fusion of one or more true groups
# by comparison to a reference map
# report all fusions whether between homeologs or not
#

import sys

refmap = sys.argv[1]
mincount = int(sys.argv[2])
minprop = float(sys.argv[3])
inpfiles = sys.argv[4:]

#load reference map
reflg = {} #which lg is marker in in the reference map
refsize = {}
f = open(refmap)
for line in f:
    tok = line.strip().split(',')
    uid = tok[0]
    lg = tok[1]
    reflg[uid] = lg
    if not lg in refsize: refsize[lg] = 0
    refsize[lg] += 1 #count markers in this ref lg
f.close()

#collect marker stats from the provisional (new) map
newlg = {}
newsize = {}
for fname in inpfiles:
    lg = fname.replace('.loc','')
    newlg[lg] = {}
    newsize[lg] = 0
    
    f = open(fname)
    for line in f:
        uid = line.strip().split()[0]
        newsize[lg] += 1 #how many markers in this new lg
        if not uid in reflg: continue
        x = reflg[uid]
        if not x in newlg[lg]: newlg[lg][x] = 0
        newlg[lg][x] += 1 #how many markers from this reflg in this new lg
    f.close()
    
for lg in newlg:
    d = newlg[lg]
    d = {x:d[x] for x in d if d[x] >= mincount} #must be a minimum count
    d = {x:d[x] for x in d if d[x] >= minprop * refsize[x]} #minimum proportion of the ref lg
    if len(d) < 2: continue
    
    print lg, ' '.join(d.iterkeys())

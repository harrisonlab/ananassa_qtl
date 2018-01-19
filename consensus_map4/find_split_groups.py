#!/usr/bin/python

#
# work out which true LGs are split into two or more provisional LGs
# by comparison to a reference map
# report all potential splits
#

import sys

refmap = sys.argv[1]
inpfiles = sys.argv[2:]

#load reference map
uid2lg = {} #which lg is marker in in the reference map
lgsize = {}
f = open(refmap)
for line in f:
    tok = line.strip().split(',')
    uid = tok[0]
    lg = tok[1]
    uid2lg[uid] = lg
f.close()

lg2frag = {}

#for each fragment (new lg) count how many hits to true lgs it has
for fname in inpfiles:
    frag = fname.replace('.loc','')
    fraghits = {}
    
    f = open(fname)
    for line in f:
        uid = line.strip().split()[0]
        if not uid in uid2lg: continue
        lg = uid2lg[uid]
        if not lg in fraghits: fraghits[lg] = 0
        fraghits[lg] += 1 #how many markers from this newlg in this reflg
    f.close()

    if len(fraghits) == 0:
        lg = 'unknown'
    else:
        cts = [[x,fraghits[x]] for x in fraghits]
        cts.sort(key=lambda x:x[1],reverse=True)
        lg = cts[0][0]
        
    if not lg in lg2frag: lg2frag[lg] = []
    lg2frag[lg].append(frag)

out = [[lg,lg2frag[lg]] for lg in lg2frag]

out.sort(key=lambda x:x[0])

for x in out: print x[0],','.join(x[1])

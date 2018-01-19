#!/usr/bin/python

'''
plot % raw-genotype-code difference of each sample wrt both parents
versus % marker pairs where neither call was missing
classify as maternal or paternal parent, progeny or nonprogeny
based on simple threshold values set by eye and entered into conf/populations
'''

import sys
import matplotlib.pyplot as plt

#figure size, inches
width = 5.0
height = 5.0
dpi = 600

#plot config
linewidth = 0.2
alpha = 0.5

#samples_vs_parents.py ${datadir}/markers/phr_and_nmh_trans.tsv ${batch}/parents ${datadir}/progeny_vs_parents.png 
datafile = sys.argv[1] #raw genotype calls
p1 = sys.argv[2]       #maternal parent
p2 = sys.argv[3]       #paternal parent
vals = [float(x) for x in sys.argv[4].split(',')] #classification thresholds
outfile = sys.argv[5]  #png file to output plot to
#classifications are output to stdout
popn = sys.argv[6]
mat_name = popn[:2]
pat_name = popn[3:]

progeny_minx = vals[0]
progeny_miny = vals[1]
p1_minx = vals[2]
p2_miny = vals[3]
sample_data = {}

f = open(datafile)
f.readline() #skip header
for line in f:
    tok = line.strip().split('\t')
    uid = tok[0]
    calls = [int(x) for x in tok[1:]]
    
    assert uid not in sample_data, "duplicated key %s"%uid
    
    sample_data[uid] = calls
f.close()

nsamples = len(sample_data)
nmarkers = len(calls)

#find full uid of parents
uid1 = p1
uid2 = p2

fig = plt.figure(figsize=(width, height)) 

ax = plt.subplot2grid((1,1), [0,0], rowspan=1)

cc = []
ss = []
xx = []
yy = []
pts = []

colour_mat   = '#ff0000'
colour_pat   = '#00ff00'
colour_pro   = '#0000ff'
colour_other = '#777777'

#compare all samples against the putative parents
for uid in sample_data:
    total = [0,0]
    diff = [0,0]
    
    for k in xrange(nmarkers):
        if sample_data[uid][k] == -1: continue
        
        if sample_data[uid1][k] != -1:
            total[0] += 1
            if sample_data[uid][k] != sample_data[uid1][k]:
                diff[0] += 1

        if sample_data[uid2][k] != -1:
            total[1] += 1
            if sample_data[uid][k] != sample_data[uid2][k]:
                diff[1] += 1
            
    if total[0] == 0 or total[1] == 0: continue
    x = 1.0 - float(diff[0]) / float(total[0])
    y = 1.0 - float(diff[1]) / float(total[1])
    
    if x >= p1_minx and y < p2_miny:
        c = colour_mat
        s = 2
        print uid,'mat',uid1,uid2
    elif x < p1_minx and y >= p2_miny:
        c = colour_pat
        s = 2
        print uid,'pat',uid1,uid2
    elif x >= progeny_minx and y >= progeny_miny:
        c = colour_pro
        s = 2
        print uid,'pro',uid1,uid2
    else:
        c = colour_other
        s = 1
        print uid,'rogue',uid1,uid2
        
    pts.append([x,y,c,s])
        
pts.sort(key=lambda x:x[3]) #make sure parents are on top
xx = [x[0] for x in pts]
yy = [x[1] for x in pts]
cc = [x[2] for x in pts]
ss = [x[3] for x in pts]

minx=min(xx)
maxx=max(xx)
miny=min(yy)
maxy=max(yy)

ax.scatter(xx,yy,c=cc,s=ss,alpha=alpha,linewidth=linewidth)
ax.vlines([progeny_minx,p1_minx],miny,maxy,alpha=0.5,linewidth=0.3)
ax.hlines([progeny_miny,p2_miny],minx,maxx,alpha=0.5,linewidth=0.3)
ax.set_xlabel(uid1 + '(%s)'%mat_name)
ax.set_ylabel(uid2 + '(%s)'%pat_name)
ax.set_title(popn)
    
plt.savefig(outfile, dpi=dpi) #, bbox_inches='tight'

#!/usr/bin/python

#
# per population:
# calc sample-vs-sample %difference
# for a random subset of markers
# plot %difference versus %missing
# collect all relicates of each sample together into a single merged sample
# output maternal and paternal samples first and second
# followed by alphabetically sorted samples
#

import sys
import os
import random
import matplotlib.pyplot as plt

datafile = sys.argv[1]    #phr_and_nmh_trans.tsv
infofile = sys.argv[2]    #file detailing which samples are male/female parent
outdir   = sys.argv[3]    #output directory
pngfile  = sys.argv[4]    #png to save plot to

mergefile = outdir + '/merged_samples'
outfile = outdir + '/mhr_phr_nmh_trans_merged.tsv'

#figure size, inches
width = 5.0
height = 5.0
dpi = 600

#choose the same marker subset each time to assist debugging
random.seed(1234)

threshold = float(os.environ['MERGE_THRESHOLD'])

#sample this many markers
nsubset = int(os.environ['MERGE_SUBSET'])

#scatter plot config
linewidth = 0.2
s = 2.0
alpha = 0.5

#load lists of maternal and paternal samples
f = open(infofile)
mat = {uid:True for uid in f.readline().strip().split()}
pat = {uid:True for uid in f.readline().strip().split()}
f.close()

data = {} #all samples

f = open(datafile)
header = f.readline()
for i,line in enumerate(f):
    tok = line.strip().split('\t')
    uid = tok[0]
    calls = [int(x) for x in tok[1:]]
    data[uid] = calls
f.close()

nsamples = len(data)
nmarkers = len(calls)

#random subset of markers
markers = range(nmarkers)
random.shuffle(markers)

fig = plt.figure(figsize=(width, height)) 

ax = plt.subplot2grid((1,1), [0,0], rowspan=1)

cc = []
xx = []
yy = []
pts = []

sample_uids = data.keys()

orig2merg = {}
merglist = {}

#all vs all
for i,uid_i in enumerate(sample_uids):
    sample_i = data[uid_i]
    
    for uid_j in sample_uids[i+1:]:
        sample_j = data[uid_j]
        
        total = 0
        diff = 0
    
        for k in markers[:nsubset]:
            #ignore missing calls
            if sample_i[k] == -1 or sample_j[k] == -1: continue
            
            total += 1
            if sample_i[k] != sample_j[k]: diff += 1
            
        x = float(total) / float(nsubset) #usable data
        y = float(diff) / float(total)    #difference
        
        if y < threshold:
            #print uid_i,uid_j #report as duplicate
            if uid_i in orig2merg and uid_j in orig2merg:
                #check both belong to same merged sample
                assert orig2merg[uid_i] == orig2merg[uid_j]
            elif uid_i in orig2merg:
                #add uid_j to existing merged sample
                mergid = orig2merg[uid_i]
                orig2merg[uid_j] = mergid
                merglist[mergid].append(uid_j)
            else:
                #this case should never happen assuming transitivity
                assert uid_j not in orig2merg

                #neither uid merged before, create new merged sample id
                mergid = '%s_merged'%uid_i
                orig2merg[uid_i] = mergid
                orig2merg[uid_j] = mergid
                merglist[mergid] = [uid_i,uid_j]
        
        pts.append([x,y])

#output scatter plot
xx = [x[0] for x in pts]
yy = [x[1] for x in pts]
ax.hlines(threshold, min(xx), max(xx), alpha=alpha)
ax.scatter(xx,yy,s=s,alpha=alpha,linewidth=linewidth)
ax.set_xlabel(datafile)
plt.savefig(pngfile, dpi=dpi)

#output list of merged samples
fout = open(mergefile,'wb')
for uid in merglist: fout.write(uid +' ' + ' '.join(merglist[uid]) + '\n')
fout.close()

#perform merging
for mergid in merglist:      #for each merged sample
    mlist = merglist[mergid]
    data[mergid] = []
    nreps = len(mlist)

    #check that merging agrees with exising mat and pat designations
    for x in mlist:
        assert (x in mat) == (mlist[0] in mat)
        assert (x in pat) == (mlist[0] in pat)

    for i in xrange(nmarkers):                #for each marker
        vals = [data[uid][i] for uid in mlist]
        vals = [x for x in vals if x != -1]

        if len(vals) == 0:
            #all calls missing
            data[mergid].append(-1)
            continue

        cts = [[x,vals.count(x)] for x in [0,1,2]]
        cts.sort(key=lambda x:x[1],reverse=True)

        if cts[0][1] > cts[1][1]:
            #assign call with unique highest count
            data[mergid].append(cts[0][0])
            continue

        #tied counts, set to missing
        data[mergid].append(-1)

#find the maternal and paternal samples, which might be merged, and might be the same
#drop replicated samples, retain only the merged sample
mat_uid = None
pat_uid = None
order = []
for uid in data:
    if uid in orig2merg: continue #skip replicated samples, only output merged version

    if uid.endswith("_merged"): tmp_uid = uid[:-7]
    else:                       tmp_uid = uid
    
    if tmp_uid in mat: mat_uid = uid #note, same uid for mat and pat in the F2 popn!
    if tmp_uid in pat: pat_uid = uid
    if not tmp_uid in mat and not tmp_uid in pat: order.append(uid)
    
assert mat != None and pat != None
order.sort()

#output in order: mat, pat, progeny by alphabetical order
fout = open(outfile,'wb')
fout.write(header)
fout.write(mat_uid + '\t' + '\t'.join([str(x) for x in data[mat_uid]]) + '\n')
fout.write(pat_uid + '\t' + '\t'.join([str(x) for x in data[pat_uid]]) + '\n')
for uid in order: fout.write(uid + '\t' + '\t'.join([str(x) for x in data[uid]]) + '\n')
fout.close()

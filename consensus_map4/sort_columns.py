#!/usr/bin/python

#
# sort columns in a loc file into a predefined order
#

import sys

locfile = sys.argv[1]
currorderfile = sys.argv[2]
neworderfile = sys.argv[3]
outfile = sys.argv[4]

#load current ordering information
currorder_dict = {}
currorder_list = []
f = open(currorderfile)
for i,line in enumerate(f):
    uid = line.strip()
    currorder_dict[uid] = i
    currorder_list.append(uid)
f.close()

#load new ordering information
neworder_list = []
neworder_dict = {}
f = open(neworderfile)
for line in f:
    uid = line.strip()
    assert uid in currorder_dict
    #if uid not in currorder_dict: continue #filter out any samples not present in the current set of columns
    neworder_dict[uid] = len(neworder_list)
    neworder_list.append(uid)
f.close()

#fout = open('popn_list_final','wb')
#for x in neworder_list: fout.write(x + '\n')
#fout.close()

n_samples = len(currorder_list)
conv = [None] * n_samples

#prepare the mapping from current column number to new column number
for i,uid in enumerate(currorder_list):
    conv[i] = neworder_dict[uid]

#load and output in reordered form, the locfile
newcalls = ['xx'] * n_samples
fout = open(outfile,'wb')
f = open(locfile)
for line in f:
    tok = line.strip().split()
    #uid,mtype,phase = tok[:3]
    calls = tok[3:]
    for i,x in enumerate(calls): newcalls[conv[i]] = x
    fout.write(' '.join(tok[:3]) + ' ' + ' '.join(newcalls) + '\n')
f.close()
fout.close()

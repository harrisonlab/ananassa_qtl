#!/usr/bin/python

'''
for each linkage group in a set of directories
pick the one with the shortest map distance as the best
assemble the best into an output directory
'''

import sys
import os
import glob

outdir = sys.argv[1]
inplist = sys.argv[2:]

#check all input directories contain the same set of lg fragments
fraglist = None
for inpdir in inplist:
    if fraglist == None:
        fraglist = [x.split('/')[-1] for x in glob.glob(inpdir+'/*.map')]
        fraglist.sort()
    else:
        tmplist = [x.split('/')[-1] for x in glob.glob(inpdir+'/*.map')]
        tmplist.sort()
        assert fraglist == tmplist

for frag in fraglist:
    best_score = 9e99
    best_dir = None
    
    for inpdir in inplist:
        f = open(inpdir + '/' + frag)
        data = []
        for line in f:
            if line.startswith('group'): continue
            tok = [x.strip() for x in line.strip().split()]
            data.append(tok)
        f.close()

        mat = [float(row[1]) for row in data if row[1] != 'NA']
        mat_max = max(mat) if len(mat) > 0 else 0.0
        
        pat = [float(row[2]) for row in data if row[2] != 'NA']
        pat_max = max(pat) if len(pat) > 0 else 0.0

        score = mat_max + pat_max
        if score < best_score:
            best_score = score
            best_dir = inpdir
            
    #print frag,best_dir,best_score
    os.system('cp %s %s'%(best_dir+'/'+frag,outdir))
    locfile = frag.replace('.map','.loc')
    os.system('cp %s %s'%(best_dir+'/'+locfile,outdir))
    

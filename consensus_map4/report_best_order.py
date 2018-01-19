#!/usr/bin/python

'''
report the sum of maternal and paternal maps for each of a set of candidate orderings
'''

import sys
import os
import glob

inplist = sys.argv[1:]

#check all input directories contain the same set of lg fragments
for inpdir in inplist:
    fraglist = [x.split('/')[-1] for x in glob.glob(inpdir+'/*.map')]

for frag in fraglist:
    print
    print frag
    best_score = 9e99
    best_dir = None
    
    for inpdir in inplist:
        try:
            f = open(inpdir + '/' + frag)
        except:
            continue
            
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
        print frag,inpdir,score
        
        if score < best_score:
            best_score = score
            best_dir = inpdir
            
    print '==>',frag,best_dir,best_score
    

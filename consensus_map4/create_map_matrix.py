#!/usr/bin/python

'''
combine the individual maps into a single matrix file suitable for loading into R
also output a mask matrix with all markers selected for prediction in all maps
'''

import os
import glob
import sys

poplist = \
[
    'RGxHA',
    'EMxFE',
    'FLxCH',
    'BSxEL',
    'HOxKO',
    'CAxDO',
    'CAxCF',
    'DAxMO',
    'P1xP2',
    'RExRE',
]

inpsubdir = sys.argv[1]
outdir = sys.argv[2]
min_markers = int(sys.argv[3]) #min markers per fragment to be included
min_map_length = float(sys.argv[4]) #min map distance of fragment to be included

lglist = []
for i in xrange(1,8):
    for x in 'ABCD':
        lglist.append(str(i)+x)
        
for lg in lglist:
    markers = {}
    all_frags = {}
    for pop in poplist:
        inpdir = 'popn_'+pop+'/map/'+inpsubdir
        globstr = inpdir+'/'+lg+'*.map'
        fraglist = glob.glob(globstr)
        
        for frag in fraglist:
            fragname = frag.split('/')[-1].replace('.map','') #eg 7A.2
            if len(fragname.split('.')) == 1:
                fragname += '.1' #eg 7A ==> 7A.1
            fragname = pop+'_'+fragname
            
            f = open(frag)
            count = 0
            for line in f:
                if not line.startswith('group'):
                    count += 1
            f.close()

            print frag,fragname,count
            if count < min_markers: continue

            f = open(frag)
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
            
            if mat_max >= min_map_length:
                #include maternal map positions
                for row in data:
                    if row[1] == 'NA': continue
                    uid = row[0]
                    pos = float(row[1])
                    if not uid in markers: markers[uid] = {}
                    markers[uid]['m'+fragname] = pos
                    
                    all_frags['m'+fragname] = None

            if pat_max >= min_map_length:
                #include paternal map positions
                for row in data:
                    if row[2] == 'NA': continue
                    uid = row[0]
                    pos = float(row[2])
                    if not uid in markers: markers[uid] = {}
                    markers[uid]['p'+fragname] = pos
                    all_frags['p'+fragname] = None
    
    outfile = outdir + '/' + lg + '.csv'
    fout = open(outfile,'wb')
    
    ordered_frags = all_frags.keys()
    ordered_frags.sort() #ensure fragments that need merging appear adjacent in the column ordering
    for frag in ordered_frags: fout.write(','+frag)
    fout.write('\n')
    
    for uid in markers:
        fout.write(uid)
        for frag in ordered_frags:
            if frag in markers[uid]:
                fout.write(',%.3f'%markers[uid][frag])
            else:
                fout.write(',')
        fout.write('\n')
    fout.close()

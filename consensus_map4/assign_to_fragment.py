#!/usr/bin/python

'''
for maps split into multiple LG fragments assign missing markers to only one fragment
and mask them out for all others so that each map uses only the nearest fragment to
predict missing marker positions
'''

import os
import glob
import sys
import numpy as np
outdir = sys.argv[1]

csvlist = glob.glob(outdir+'/*.csv')

for csvfile in csvlist:
    lg = csvfile.split('/')[-1].replace('.csv','')
    
    #load consensus map position for all markers
    mapfile = csvfile.replace('.csv','.out')
    conpos = {}
    f = open(mapfile)
    f.readline() #skip header
    for line in f:
        tok = line.strip().split(',')
        uid = tok[0].strip('"')
        pos = float(tok[1])
        conpos[uid] = pos 
    f.close()
    
    #group lg fragments by map
    f = open(csvfile)
    header = f.readline().strip().split(',')
    frags = {}
    for x in header[1:]:
        mp = x.split('_')[0]
        if not mp in frags: frags[mp] = []
        frags[mp].append(x)
        
    #retain only lgs with more than one fragment
    frags = {x:frags[x] for x in frags if len(frags[x]) > 1}
        
    #find which fragments are part of split LGs
    splitlgs = {}
    for mp in frags:
        for x in frags[mp]:
            splitlgs[x] = []
        
    #determine which markers are present in which fragments
    marker = {}
    for line in f:
        tok = line.strip().split(',')
        uid = tok[0]
        posn = tok[1:]
        for i,x in enumerate(header[1:]):
            mp = x.split('_')[0]
            if not uid in marker: marker[uid] = {}
            if posn[i] != '': marker[uid][x] = None #we only care if marker is present or absent
    
    f.close()
    
    #find min, median and max consensus marker position per fragment
    stats = {}
    for x in splitlgs:
        for uid in marker:
            if x in marker[uid]:
                splitlgs[x].append(conpos[uid])
        stats[x] = [min(splitlgs[x]),np.median(splitlgs[x]),max(splitlgs[x])]
        
    #assign all markers to only one fragment per map
    assignment = {}
    for uid in marker:
        for mp in frags:
            if not uid in assignment: assignment[uid] = {}
            assignment[uid][mp] = None
            
            #marker should be either absent from all fragments
            #or present in only one
            
            #is marker absent from all fragments?
            location = None
            for x in frags[mp]:
                if x in marker[uid]:
                    assert location == None #check not present in more than one fragment
                    location = x
                
            #assign marker to fragment that contains it
            if location != None:
                assignment[uid][mp] = location
                continue
                
            #work out which fragment is closest
            flist = []
            mposn = conpos[uid]
            
            #is the marker contained within one fragment only?
            for x in frags[mp]:
                if mposn >= stats[x][0] and mposn <= stats[x][2]:
                    flist.append(x)
            
            if len(flist) == 1:
                assignment[uid][mp] = flist[0]
                continue

            #assign to fragment with closest median position
            flist = []
            for x in frags[mp]: flist.append([x,abs(mposn-stats[x][1])])
            flist.sort(key=lambda y:y[1])
            assignment[uid][mp] = flist[0][0]
        
    #create mask file
    f = open(csvfile)
    header = f.readline()

    maskfile = csvfile+'.mask'
    fout = open(maskfile,'wb')
    fout.write(header)
    
    fraglist = header.strip().split(',')[1:]
    maplist = [x.split('_')[0] for x in fraglist]
    
    for line in f:
        tok = line.strip().split(',')
        uid = tok[0]
        fout.write(uid)
        
        for i in xrange(len(fraglist)):
            mp = maplist[i]
            fg = fraglist[i]
            
            if uid not in assignment or mp not in assignment[uid]:
                fout.write(',FALSE') #do not mask
                continue
                
            if assignment[uid][mp] == fg:
                fout.write(',FALSE') # do not mask
            else:
                fout.write(',TRUE') # mask
            
        fout.write('\n')

    fout.close()
    f.close()

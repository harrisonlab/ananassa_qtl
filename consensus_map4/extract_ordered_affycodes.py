#!/usr/bin/python
#Crosslink Copyright (C) 2016 NIAB EMR see included NOTICE file for details

#
# convert from two columns per sample to one column per sample
# preserver the ordering of the allele codes eg A,B, -> AB,

import sys

inpfile = sys.argv[1]

fout = sys.stdout

#marker,cM,type,phase,mat,mat,pat,pat,RGxHA1001-mat,RGxHA1001-pat,RGxHA002-mat...
f = open(inpfile)
header = f.readline().strip().split(',')
fout.write('marker,type,phase,maternal,paternal')

for i,x in enumerate(header[8:]):
    if i%2 == 1: continue
    fout.write(','+x.split('-')[0].strip())

fout.write('\n')

for line in f:
    tok = line.strip().split(',')
    marker = tok[0]
    pos = tok[1]
    mtype = tok[2]
    phase = tok[3]
    mat0 = tok[4]
    mat1 = tok[5]
    pat0 = tok[6]
    pat1 = tok[7]

    calls = tok[8:]

    fout.write(','.join([marker,mtype,phase,mat0+mat1,pat0+pat1]))

    for i,x in enumerate(calls):
        if i%2 == 0: fout.write(',')
        fout.write(x)

    fout.write('\n')
f.close()


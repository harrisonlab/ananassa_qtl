#!/usr/bin/python

#
# output haplotype information
#

import sys

locfile = sys.argv[1]
mapfile = sys.argv[2]
parfile = sys.argv[3] #original parental genotype calls
profile = sys.argv[4] #progeny names
outbase = sys.argv[5]

#load parental genotype calls
parcall = {}
f = open(parfile)
for line in f:
    tok = line.strip().split('\t')
    uid = tok[0]
    mat = tok[1]
    pat = tok[2]
    parcall[uid] = [mat,pat]
f.close()

#load progeny names
f = open(profile)
prog = [line.strip() for line in f]
f.close()

#load map positions
posn = {}
f = open(mapfile)
for line in f:
    uid,lg,pos = line.strip().split(',')
    posn[uid] = float(pos)
f.close()

#process genotype information
conv = {'2':'0','0':'2'}
conv2 = {'0':'AA', '1':'AB', '2':'BB'}
conv3 = {'A':'B','B':'A'}

marker = {}
f = open(locfile)
for line in f:
    tok = line.strip().split()
    uid,mtype,phase = tok[0:3]
    
    call = tok[3:]
    assert len(call) == len(prog)
    
    mat = parcall[uid][0]
    pat = parcall[uid][1]
    
    #assuming these markers have been filtered out
    #apply type correction if required
    #if mtype == '<nnxnp>' and mat == '1':   ##eg 12  => 01
    #    mat = conv[pat]
    #    pat = '1'
    #elif mtype == '<lmxll>' and pat == '1':
    #    pat = conv[mat]
    #    mat = '1'
    
    #convert from numerical to explicit AB codes
    if mat == '1' and pat == '0': #lmxll
        mat = 'AB'
        pat = 'AA'
    elif mat == '1' and pat == '2':
        mat = 'BA'
        pat = 'BB'
    elif mat == '0' and pat == '1':#nnxnp
        mat = 'AA'
        pat = 'AB'
    elif mat == '2' and pat == '1':
        mat = 'BB'
        pat = 'BA'
    else: #hkxhk
        mat = 'AB'
        pat = 'AB'
    
    #apply phasing
    if phase[1] == '1': mat = mat[::-1]
    if phase[2] == '1': pat = pat[::-1]
    
    #convert calls into AB code
    newcall = []
    hapcode = []
    for x in call:
        assert not '-' in x
        if x == 'lm':
            #heterozygote
            m = mat[1]
            p = pat[0]
            hm = 'X'
            hp = '-'
        elif x == 'll':
            #homozygote
            m = mat[0]
            p = pat[0]
            hm = 'O'
            hp = '-'
        elif x == 'np':
            #heterozygote
            m = mat[0]
            p = pat[1]
            hm = '-'
            hp = 'X'
        elif x == 'nn':
            #homozygote
            m = mat[0]
            p = pat[0]
            hm = '-'
            hp = 'O'
        else:
            if x[0] == 'h':
                m = 'A'
            else:
                m = 'B'
            if x[1] == 'h':
                p = 'A'
            else:
                p = 'B'
            hm = 'O' if mat[0] == m else 'X'
            hp = 'O' if pat[0] == p else 'X'
            hm = '-'
            hp = '-'
            
        newcall.append(m+p)
        hapcode += [hm,hp]
                                                    #6
    marker[uid] = [uid,posn[uid],mtype,phase,mat,pat,call,newcall,hapcode]
f.close()

data = []
for uid in marker: data.append(marker[uid])
data.sort(key=lambda row:row[1])

#data = [row for row in data if row[1] == '<lmxll>']

#output crosslink genotypes
fout = open(outbase + '_genotypes_crosslink.csv','wb')
fout.write('marker,cM,type,phase,mat,pat,' + ','.join(prog) + '\n')
for row in data:
    mtype = row[2]
    tmprow = row[:6]
    tmprow[4] = mtype[1:3]
    tmprow[5] = mtype[4:6]
    fout.write(','.join([str(x) for x in tmprow]) + ',')
    fout.write(','.join(row[6]) + '\n')
fout.close()

#output affy genotypes
fout = open(outbase + '_genotypes_affy.csv','wb')
fout.write('marker,cM,type,phase,mat,pat,' + ','.join(prog) + '\n')
for row in data:
    mtype = row[2]
    tmprow = row[:6]
    tmprow[4] = row[4]
    tmprow[5] = row[5]
    fout.write(','.join([str(x) for x in tmprow]) + ',')
    fout.write(','.join(row[6]) + '\n')
fout.close()

#output haplotypes
fout = open(outbase + '_haplotypes.csv','wb')
fout.write(',,,mat,pat\n')
for row in data:
    fout.write(','.join([str(x) for x in row[:4]]))
    fout.write(',' + ','.join(row[6]) + '\n')
fout.close()

#output haplocodes
fout = open(outbase + '_haplocodes.csv','wb')
fout.write(',,,mat,pat\n')
for row in data:
    fout.write(','.join([str(x) for x in row[:7]]))
    fout.write(',' + ','.join(row[7]) + '\n')
fout.close()

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

#load parental genotype calls as affy codes
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

marker = {}
f = open(locfile)
for line in f:
    tok = line.strip().split()
    uid,mtype,phase = tok[0:3]

    call = tok[3:]
    assert len(call) == len(prog)

    mat = parcall[uid][0]
    pat = parcall[uid][1]
    matpat = mat+pat

    #generate haplotype code (identify the two parental haplotypes as 0 and 1)
    hapcode = []

    if mtype == '<lmxll>':
        if phase[1] == '0': conv = {'lm':'1','ll':'0'}
        else:               conv = {'lm':'0','ll':'1'}

        for x in call: hapcode += [conv[x],'-']
    elif mtype == '<nnxnp>':
        if phase[2] == '0': conv = {'np':'1','nn':'0'}
        else:               conv = {'np':'0','nn':'1'}

        for x in call: hapcode += ['-',conv[x]]

    elif mtype == '<hkxhk>':
        if phase[1] == '0': conv1 = {'hh':'0','hk':'0','kh':'1','kk':'1'}
        else:               conv1 = {'hh':'1','hk':'1','kh':'0','kk':'0'}
        if phase[2] == '0': conv2 = {'hh':'0','hk':'1','kh':'0','kk':'1'}
        else:               conv2 = {'hh':'1','hk':'0','kh':'1','kk':'0'}

        for x in call: hapcode += [conv1[x],conv2[x]]

    #work out parental haplotypes
    if mtype == '<lmxll>':
        if matpat == '10':
            mat = 'BA' if phase[1] == '0' else 'AB'
            pat = 'AA'
        else:
            assert matpat == '12'
            mat = 'AB' if phase[1] == '0' else 'BA'
            pat = 'BB'

    elif mtype == '<nnxnp>':
        if matpat == '01':
            pat = 'BA' if phase[2] == '0' else 'AB'
            mat = 'AA'
        else:
            assert matpat == '21'
            pat = 'AB' if phase[2] == '0' else 'BA'
            mat = 'BB'

    else: #hkxhk
        mat = 'AB' if phase[1] == '1' else 'BA'
        pat = 'AB' if phase[2] == '1' else 'BA'

    #generate haplotypes as affy A/B allele codes
    hap = []
    conv = {'A':'B', 'B':'A', '-':'-'} #switch A <=> B
    if mtype == '<lmxll>':
        for x in call:
            if x == 'lm': hap += [conv[pat[0]],pat[0]] #het
            else:         hap += [pat[0],pat[0]]       #hom

    elif mtype == '<nnxnp>':
        for x in call:
            if x == 'np': hap += [mat[0],conv[mat[0]]] #het
            else:         hap += [mat[0],mat[0]]       #hom

    elif mtype == '<hkxhk>':
        for x in call:
            if x == 'hh':   hap += ['A','A']
            elif x == 'kk': hap += ['B','B']
            elif x == 'kh': hap += ['B','A']
            else:           hap += ['A','B']
                                             #4   5   6    7   8
    marker[uid] = [uid,posn[uid],mtype,phase,mat,pat,call,hap,hapcode]
f.close()

data = []
for uid in marker: data.append(marker[uid])
data.sort(key=lambda row:row[1])

#data = [row for row in data if row[3] != '<lmxll>']

#output genotypes
fout = open(outbase + '_genotypes.csv','wb')
fout.write('marker,cM,type,phase,mat,pat,' + ','.join(prog) + '\n')
for row in data:
    mtype = row[2]
    tmprow = row[:6]
    tmprow[4] = mtype[1:3]
    tmprow[5] = mtype[4:6]
    fout.write(','.join([str(x) for x in tmprow]) + ',')
    fout.write(','.join(row[6]) + '\n')
fout.close()

#output haplotype
fout = open(outbase + '_haplotypes.csv','wb')
fout.write('marker,cM,type,phase,mat,mat,pat,pat,' + ','.join([x+'-mat,'+x+'-pat' for x in prog]) + '\n')
for row in data:
    fout.write(','.join([str(x) for x in row[:4]]))
    fout.write(',' + ','.join([row[4][0],row[4][1],row[5][0],row[5][1]]))
    fout.write(',' + ','.join(row[7]) + '\n')
fout.close()

#output haplocodes
fout = open(outbase + '_haplocodes.csv','wb')
fout.write('marker,cM,type,phase,mat,mat,pat,pat,' + ','.join([x+'-mat,'+x+'-pat' for x in prog]) + '\n')
for row in data:
    fout.write(','.join([str(x) for x in row[:4]]))
    fout.write(',' + ','.join([row[4][0],row[4][1],row[5][0],row[5][1]]))
    fout.write(',' + ','.join(row[8]) + '\n')
fout.close()

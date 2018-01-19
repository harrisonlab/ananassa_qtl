#!/usr/bin/python

'''
recode markers from Affymetrix 0 (AA), 1 (AB), 2 (BB) and -1 (NoCall)
into joinmap compatible coding
output a dummy phase (always zero)
output xx where the allele is expected to be impossible given the parental genotypes
modify the marker name to include the maternal and paternal affy code
'''

import sys

fname = sys.argv[1]
popn = sys.argv[2]

f = open(fname)
f.readline() #skip header

#count markers by type
ct = {}
type_list = ['MissingParental','Monomorphic','<lmxll>','<nnxnp>','<hkxhk>']
for x in type_list: ct[x] = 0

for line in f:
    tok = line.strip().split('\t')
    
    name = tok[0]
    
    calls = [int(x) for x in tok[1:]]

    mat = calls[0]
    pat = calls[1]
    calls = calls[2:]
    
    if mat == 1 and pat == 1:
        #hkxhk
        mtype = '<hkxhk>'
        conv = {0:'hh',1:'hk',2:'kk',-1:'--'}
        aux = '11'
        
    elif mat == 1 and pat == 0:
        #lmxll
        mtype = '<lmxll>'
        conv = {0:'ll',1:'lm',2:'xx',-1:'--'}
        aux = '10'

    elif mat == 1 and pat == 2:
        #lmxll
        mtype = '<lmxll>'
        conv = {0:'xx',1:'lm',2:'ll',-1:'--'}
        aux = '12'

    elif mat == 0 and pat == 1:
        #nnxnp
        mtype = '<nnxnp>'
        conv = {0:'nn',1:'np',2:'xx',-1:'--'}
        aux = '01'

    elif mat == 2 and pat == 1:
        #nnxnp
        mtype = '<nnxnp>'
        conv = {0:'xx',1:'np',2:'nn',-1:'--'}
        aux = '21'
    elif mat == -1 and pat == -1:
        mtype = 'MissingParental'
    elif mat == -1:
        mtype = 'MissingParental'
    elif pat == -1:
        mtype = 'MissingParental'
    elif (mat == 0 and pat == 0) or (mat == 2 and pat == 2):
        mtype = 'Monomorphic'
    elif (mat == 0 and pat == 2) or (mat == 2 and pat == 0):
        mtype = 'Monomorphic'
    else:
        #other
        assert False
        
    ct[mtype] += 1
    
    #drop unusable markers
    if mtype not in ['<lmxll>','<nnxnp>','<hkxhk>']: continue
    
    calls = [conv[x] for x in calls]
    
    if mtype == '<lmxll>': phase = '{0-}'
    if mtype == '<nnxnp>': phase = '{-0}'
    if mtype == '<hkxhk>': phase = '{00}'
    
    print name[:3] + aux + name[3:], mtype, phase, ' '.join(calls)

#write out marker type counts

keys = [x for x in ct if x.startswith('<')]
ct['Informative'] = sum([ct[x] for x in keys])
for x in keys: del ct[x]

fout = open(fname + '.out','wb')
for mtype in ct: fout.write(popn + ' ' + mtype + ' ' + str(ct[mtype]) + '\n')
fout.close()

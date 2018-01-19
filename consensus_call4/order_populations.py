#!/usr/bin/python

#
# extract data for one population from the raw genotype calls
#

import sys

classfile = sys.argv[1]
callfile = sys.argv[2]
outdata = sys.argv[3]
outinfo = sys.argv[4]

mat = {}
pat = {}
pro = {}

f = open(classfile)
for line in f:
    tok = line.strip().split()
    uid = tok[0]
    typ = tok[1]

    if typ == 'mat':     mat[uid] = None
    elif typ == 'pat':   pat[uid] = None
    elif typ == 'pro':   pro[uid] = None
    elif typ == 'rogue': pass
    else:                assert False, 'unknown type %s'%typ
f.close()

data = {}

f = open(callfile)
header = f.readline()
for line in f:
    tok = line.strip().split()
    uid = tok[0]
    calls = tok[1:]
    data[uid] = calls
f.close()

fout = open(outinfo,'wb')
fout.write(' '.join(mat.iterkeys()) + '\n') #maternal samples
fout.write(' '.join(pat.iterkeys()) + '\n') #paternal samples
#fout.write(' '.join(pro.iterkeys()) + '\n') #progeny samples
fout.close()

fout = open(outdata,'wb')
fout.write(header)
for uid in mat.iterkeys(): fout.write(uid + '\t' + '\t'.join(data[uid]) + '\n')
for uid in pat.iterkeys(): fout.write(uid + '\t' + '\t'.join(data[uid]) + '\n')
for uid in pro.iterkeys(): fout.write(uid + '\t' + '\t'.join(data[uid]) + '\n')
fout.close()

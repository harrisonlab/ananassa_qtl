#!/usr/bin/python

#
# process CAxDO population:
# anything with Hibrid in the name is the F1 hybrid
# exclude maternal and paternal samples
# retain other progeny as legitimate
#

import sys

classfile = sys.argv[1]
callfile = sys.argv[2]
outdata = sys.argv[3]
outinfo = sys.argv[4]

f1 = {}
pro = {}

f = open(classfile)
for line in f:
    tok = line.strip().split()
    uid = tok[0]
    typ = tok[1]

    if 'Hibrid' in uid:
        f1[uid] = None
    elif typ == 'pro':
        pro[uid] = None
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
fout.write(' '.join(f1.iterkeys()) + '\n') #maternal samples
fout.write(' '.join(f1.iterkeys()) + '\n') #paternal samples
#fout.write(' '.join(pro.iterkeys()) + '\n') #progeny samples
fout.close()

fout = open(outdata,'wb')
fout.write(header)
for uid in f1.iterkeys(): fout.write(uid + '\t' + '\t'.join(data[uid]) + '\n')
for uid in pro.iterkeys(): fout.write(uid + '\t' + '\t'.join(data[uid]) + '\n')
fout.close()

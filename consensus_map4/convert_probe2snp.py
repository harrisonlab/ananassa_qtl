#!/usr/bin/python
#Crosslink Copyright (C) 2016 NIAB EMR see included NOTICE file for details
#
# convert probeset ids to snp ids in a loc file
#

import sys

inp = sys.argv[1]
p2sname = sys.argv[2]
out = sys.argv[3]

p2s = {}
f = open(p2sname)
f.readline()
for line in f:
    tok = line.strip().split()
    pid = tok[0][3:]
    sid = tok[1]
    p2s[pid] = sid
f.close()

f = open(inp)
fout = open(out,'wb')
for line in f:
    tok = line.strip().split(' ')
    
    assert tok[0].startswith('NMH') or tok[0].startswith('PHR')
    pid = tok[0][6:]
    
    assert pid in p2s, pid
    tok[0] = p2s[pid]
            
    fout.write(' '.join(tok) + '\n')
f.close()
fout.close()

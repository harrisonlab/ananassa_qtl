#!/usr/bin/python

#
# convert probeset ids to snp ids in a csv
# usage: cat inpfile | probe2snp.py > outfile
#

import sys
import os

p2s = {}

f1 = '~/octoploid_mapping/axiom_chip_info/IStraw90.r1.ps2snp_map.ps.fixed'
f2 = '~/rjv_mnt/cluster/octoploid_mapping/axiom_chip_info/IStraw90.r1.ps2snp_map.ps.fixed'

if os.path.isfile(f1):
    f = open(f1)
else:
    f = open(f2)
    
f.readline() #skip header
for line in f:
    tok = line.strip().split()
    pid = tok[0]
    sid = tok[1]
    p2s[pid] = sid
f.close()

for line in std.stdin:
    if not line.startswith('AX-'):
        sys.stdout.write(line)
        continue
        
    pid = line[:11]          # eg AX-89803679
    sid = p2s[pid]
    line = sid + line[11:]
    sys.stdout.write(line)
    
f.close()

#!/usr/bin/python

#
# convert label1 label2 data list into label1 x label2 data table
#

import sys

data = {}
f = open(sys.argv[1])
for line in f:
    tok = line.strip().split()
    assert len(tok) == 3
    if not tok[0] in data: data[tok[0]] = {}
    data[tok[0]][tok[1]] = tok[2]
f.close()

for i,x in enumerate(data):
    if i == 0: print ',' + ','.join(data[x].iterkeys())
    print x + ',' + ','.join(data[x].itervalues())

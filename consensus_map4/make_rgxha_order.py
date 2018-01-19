#!/usr/bin/python

#
# order the sample names into "population" order
#

import sys

inp = sys.argv[1]

f = open(inp)
data = [line.strip() for line in f]
f.close()

data.sort(key=lambda x:x[-3:])

for x in data: print x

#!/usr/bin/python

#
# get just the progeny names from the tsv file
#

import sys

inp = sys.argv[1]

f = open(inp)
header = f.readline().strip().split('\t')[3:]
f.close()

for x in header: print x

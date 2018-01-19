#!/usr/bin/python

#
# merge split LGs
#


import sys
import pandas as pd

inpfile = sys.argv[1]
mapfile = sys.argv[2]

df = pd.DataFrame.from_csv(inpfile)

cols = list(df.columns.values)
cols.sort()

print inpfile
print len(cols)
print cols

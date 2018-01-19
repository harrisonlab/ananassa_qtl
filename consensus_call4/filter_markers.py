#!/usr/bin/python

'''
use chisquared to segregation filter markers at requested p value
(use pvalue of 0.0 to disable this filter)
also filter on proportion of missing data
(use pmissing of 1.0 to disable this filter)
'''

import sys
from scipy.stats import chisquare
from rjv_stats import gtest

pvalue   = float(sys.argv[1])
pmissing = float(sys.argv[2])
fname    = sys.argv[3]

f = open(fname)
f.readline() #skip header

for line in f:
    tok = line.strip().split(' ')
    
    name = tok[0]
    mtype = tok[1]
    phase = tok[2]
    calls = tok[3:]

    if mtype == '<hkxhk>':
        aa = calls.count('hh')
        ab = calls.count('hk')
        bb = calls.count('kk')
        
        obs = [aa,ab,bb]
        total = float(sum(obs))
        exp = [total*0.25,total*0.5,total*0.25]
        
    else:
        if mtype == '<lmxll>':
            aa = calls.count('ll')
            ab = calls.count('lm')
        else:
            assert(mtype == '<nnxnp>')
            aa = calls.count('nn')
            ab = calls.count('np')

        obs = [aa,ab]
        total = float(sum(obs))
        exp = [total*0.5,total*0.5]
        
    #chisq2,p2 = chisquare(obs,exp)
    
    chisq,p = gtest(obs,exp)
    
    #sys.stderr.write("%f %f     %f %f\n"%(chisq,p,chisq2,p2))
    
    missing = (float(len(calls)) - total) / float(len(calls))
    
    if p < pvalue: continue
    
    if missing >= pmissing: continue

    print line.strip()

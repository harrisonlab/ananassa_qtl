import sys
import math
from scipy.stats import chisqprob

def gtest(obs, exp, ddof=0):
    '''
    http://en.wikipedia.org/wiki/G-test
    test for goodness of fit to expected frequencies
    
    obs - observed fres
    exp - expected freqs
    ddof - delta dof

    returns
    chisquare statistic and p value
    
    based on https://gist.github.com/brentp/570896
    '''

    assert len(obs) == len(exp) 
    assert not 0.0 in exp

    n = len(obs)
    
    g = 0.0
    for i in xrange(n):
        if obs[i] == 0.0: continue #Oi * ln( Oi / Ei) == 0 if Oi == 0
        
        g += obs[i] * math.log(obs[i] / exp[i])
        
        if exp[i] < 5.0: sys.stderr.write("warning: expected value less than 5 in gtest\n")
    g *= 2.0
    
    return g, chisqprob(g, n - 1 - ddof)

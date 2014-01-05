#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file covar.m

 Python version by RDL 14 Feb 2012
 Copyright notice from covar.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""

import numpy as np
import convm 

def covar(x, p):
    """Generates a covariance  matrix
    
    Usage: X = covar(x,p)
    Generates a p x p covariance matrix for the sequence x.

    """
    x = x[:]
    m = len(x)
    x = x - np.ones(m) * np.sum(x) / m

    R = np.dot(convm.convm(x, p).conj().T, convm.convm(x, p)) / (m-1)
    
    return R

#x = x(:);
#m = length(x);
#x = x - ones(m,1)*(sum(x)/m);
#R = convm(x,p)'*convm(x,p)/(m-1);
#end;
    
def main():
    """Just a test driver, compare with Hayes pp. 572-573"""
    x = np.array([1, 3, 2])
    X = covar(x,4)
    print(X)

#    X = np.cov(x, rowvar=0)
#    print(X)
    
if __name__ == '__main__':
    main()


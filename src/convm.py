#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file convm.m

 Python version by RDL 10 Jan 2012
 Copyright notice from convm.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""

import numpy as np

def convm(x, p):
    """Generates a convolution matrix
    
    Usage: X = convm(x,p)
    Given a vector x of length N, an N+p-1 by p convolution matrix is
    generated of the following form:
              |  x(0)  0      0     ...      0    |
              |  x(1) x(0)    0     ...      0    |
              |  x(2) x(1)   x(0)   ...      0    |
         X =  |   .    .      .              .    |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |  x(N) x(N-1) x(N-2) ...  x(N-p+1) |
              |   0   x(N)   x(N-1) ...  x(N-p+2) |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |   0    0      0     ...    x(N)   |
         
    That is, x is assumed to be causal, and zero-valued after N.

    """
    N = len(x) + 2*p - 2
    xpad = np.concatenate([np.zeros(p-1), x[:], np.zeros(p-1)])
    X = np.zeros((len(x)+p-1, p))
    # Construct X column by column
    for i in xrange(p):
        X[:,i] = xpad[p-i-1:N-i]
    
    return X
    
def main():
    """Just a test driver, compare with Hayes pp. 572-573"""
    x = np.array([1, 3, 2])
    X = convm(x,4)
    print(X)
    
if __name__ == '__main__':
    main()


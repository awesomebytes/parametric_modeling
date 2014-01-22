#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file convm.m

 Python version by RDL 10 Jan 2012
 Copyright notice from convm.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""

import numpy as np

def convmtx(v, n):
    """Generates a convolution matrix
    
    Usage: X = convm(v,n)
    Given a vector v of length N, an N+n-1 by n convolution matrix is
    generated of the following form:
              |  v(0)  0      0     ...      0    |
              |  v(1) v(0)    0     ...      0    |
              |  v(2) v(1)   v(0)   ...      0    |
         X =  |   .    .      .              .    |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |  v(N) v(N-1) v(N-2) ...  v(N-n+1) |
              |   0   v(N)   v(N-1) ...  v(N-n+2) |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |   0    0      0     ...    v(N)   |
    And then it's trasposed to fit the MATLAB return value.     
    That is, v is assumed to be causal, and zero-valued after N.

    """
    N = len(v) + 2*n - 2
    xpad = np.concatenate([np.zeros(n-1), v[:], np.zeros(n-1)])
    X = np.zeros((len(v)+n-1, n))
    # Construct X column by column
    for i in xrange(n):
        X[:,i] = xpad[n-i-1:N-i]
    
    return X
    
def main():
    """Just a test"""
    h = [1,2,3,2,1]
    X = convmtx(h,7)
    print(X)
    
## MATLAB OUTPUT:
# >> h = [1 2 3 2 1];
# >> convmtx(h,7)
# 
# ans =
# 
#      1     2     3     2     1     0     0     0     0     0     0
#      0     1     2     3     2     1     0     0     0     0     0
#      0     0     1     2     3     2     1     0     0     0     0
#      0     0     0     1     2     3     2     1     0     0     0
#      0     0     0     0     1     2     3     2     1     0     0
#      0     0     0     0     0     1     2     3     2     1     0
#      0     0     0     0     0     0     1     2     3     2     1
## PYTHON OUTPUT:
# array([[ 1.,  2.,  3.,  2.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
#        [ 0.,  1.,  2.,  3.,  2.,  1.,  0.,  0.,  0.,  0.,  0.],
#        [ 0.,  0.,  1.,  2.,  3.,  2.,  1.,  0.,  0.,  0.,  0.],
#        [ 0.,  0.,  0.,  1.,  2.,  3.,  2.,  1.,  0.,  0.,  0.],
#        [ 0.,  0.,  0.,  0.,  1.,  2.,  3.,  2.,  1.,  0.,  0.],


    
if __name__ == '__main__':
    main()


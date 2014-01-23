#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Created on Jan 22 21:38 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com
This file pretends to imitate the behavior of the MATLAB function convmtx

"""

import numpy as np
from scipy.linalg import toeplitz

def convmtx(v, n):
    """From MATLAB:
    %CONVMTX Convolution matrix.
    %   CONVMTX(C,N) returns the convolution matrix for vector C.
    %   If C is a column vector and X is a column vector of length N,
    %   then CONVMTX(C,N)*X is the same as CONV(C,X).
    %   If R is a row vector and X is a row vector of length N,
    %   then X*CONVMTX(R,N) is the same as CONV(R,X).
    """
    # Local Variables: cidx, c, x_left, ridx, m, n, x_right, mv, t, v, x, r, nv
    # Function calls: convmtx, length, ones, zeros, size, toeplitz
    #%CONVMTX Convolution matrix.
    #%   CONVMTX(C,N) returns the convolution matrix for vector C.
    #%   If C is a column vector and X is a column vector of length N,
    #%   then CONVMTX(C,N)*X is the same as CONV(C,X).
    #%   If R is a row vector and X is a row vector of length N,
    #%   then X*CONVMTX(R,N) is the same as CONV(R,X).
    #%
    #%   % Example:
    #%   %   Generate a simple convolution matrix.
    #%
    #%   h = [1 2 3 2 1];
    #%   convmtx(h,7)        % Convolution matrix
    #%
    #%   See also CONV.
    #%   Author(s): L. Shure, 47-88
    #%          T. Krauss, 3-30-93, removed dependence on toeplitz
    #%   Copyright 1988-2004 The MathWorks, Inc.
    #%   $Revision: 1.6.4.3 $  $Date: 2012/10/29 19:30:54 $
    try:
        [nv, mv] = v.shape # if its vertical, shape will return 2 values, rows and cols
    except ValueError: # if its horizontal only len value will be available
        mv = len(v)
        nv = 1
    v = v.flatten(1)

    #c = np.vstack((v, np.zeros(n-1)))
    c = np.hstack((v, np.zeros((n-1))))
    r = np.zeros(n)
    m = len(c)
    x_left = r[n:0:-1] # reverse order from n to 2 in original code
    x_right = c.flatten(1)
    x = np.hstack((x_left, x_right))
    #%x = [r(n:-1:2) ; c(:)];                 % build vector of user data
    cidx = np.arange(0., (m-1.)+1).conj().T
    ridx = np.arange(n, (1.)+(-1.), -1.)

    t = np.zeros([len(cidx),len(ridx)])
    counter_cidx = 0
    for c_val in cidx:
        counter_ridx = 0
        for r_val in ridx:
            t[counter_cidx, counter_ridx] = c_val + r_val
            counter_ridx += 1
        counter_cidx += 1
    #t = cidx[:,int(np.ones(n))-1] + ridx[int(np.ones(m))-1,:] # that double loop should do this...
    #% Toeplitz subscripts

    t[:] = x[t.astype(int)-1]
    #% actual data
    #% end of toeplitz code

    if mv<nv:
        t = t.T
    
    return t
    
    
    
    
    
    

#         """Generates a convolution matrix
#     
#     Usage: X = convm(v,n)
#     Given a vector v of length N, an N+n-1 by n convolution matrix is
#     generated of the following form:
#               |  v(0)  0      0     ...      0    |
#               |  v(1) v(0)    0     ...      0    |
#               |  v(2) v(1)   v(0)   ...      0    |
#          X =  |   .    .      .              .    |
#               |   .    .      .              .    |
#               |   .    .      .              .    |
#               |  v(N) v(N-1) v(N-2) ...  v(N-n+1) |
#               |   0   v(N)   v(N-1) ...  v(N-n+2) |
#               |   .    .      .              .    |
#               |   .    .      .              .    |
#               |   0    0      0     ...    v(N)   |   
#     That is, v is assumed to be causal, and zero-valued after N.
# 
#     """
    
#     N = len(v) + 2*n - 2
#     xpad = np.concatenate([np.zeros(n-1), v[:], np.zeros(n-1)])
#     X = np.zeros((len(v)+n-1, n))
#     # Construct X column by column
#     for i in xrange(n):
#         X[:,i] = xpad[n-i-1:N-i]
#         

    #t = toeplitz(np.vstack((v, np.zeros(n-1))), np.zeros(n))
#     toep_left = np.hstack((v, np.zeros((n-1))))
#     print toep_left
#     toep_right = np.zeros(n)
#     print toep_right
    #t = toeplitz(np.array(np.vstack((np.hstack((v)), np.hstack((np.zeros((n-1.))))))), np.zeros(n))
#     t = toeplitz(np.hstack((v, np.zeros((n-1)))), np.hstack((np.zeros(n))))
#     #H = toeplitz(h / c, np.array(np.hstack((1, np.zeros(K)))))
#     try:
#         [nv, mv] = t.shape
#     except ValueError:
#         mv = len(v)
#         nv = 1
#     if mv < nv:
#         return t.T
#     return t
#    
# #     if mv < nv:
# #         return X.T
# #     return X
    
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


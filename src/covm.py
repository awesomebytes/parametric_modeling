#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file covm.m

 Python version by RDL 31 Mar 2012
 Copyright notice from covm.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""
from __future__ import print_function,division
import numpy as np
from convm import convm

def covm(x, p):
    """ Find an all-pole model using the covariance method

    USAGE:	a,err = covm(x,p) 
   
   	An all-pole of order p is found for the input sequence 
   	x using the covariance method.  The model is of the form
   		H(z) = b(0)/A(z) 
   	The coefficients of A(z) are returned in the vector
   		a=[1, a(1), ... a(p)]
   	and the modeling error is returned in err.

    """
    x = x.flatten()
    N = len(x)
    if p > N:
        print('ERROR: model order too large')
    else:
        X = convm(x, p+1)
        Xq  = X[p-1:N-1,0:p]
        xq1 = -X[p:N, 0]
        a = np.linalg.lstsq(Xq, xq1)[0]
        a = np.insert(a, 0, 1)
        err = np.dot(X[p:N,0].conj().T, X[p:N,:])
        # Added by RDL: it seems that in order to have the sample
        # autocorrelation, you have to normalize by the number of samples used,
        # which in this case is (N-p)
        err /= (N-p)
        err = np.dot(err, a)
        err = np.abs(err)
    return a, err

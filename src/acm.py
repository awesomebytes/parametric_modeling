#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file acm.m

 Python version by RDL 29 Mar 2012
 Copyright notice from acm.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""
from __future__ import print_function,division
import numpy as np
from convm import convm

def acm(x,p):
    """ Find an all-pole model using the autocorrelation method

    Usage: a,err = acm(x,p) 

	The input sequence x is modeled as the unit sample response of
	a filter having a system function of the form
		H(z) = b(0)/A(z) 
	where the coefficients of A(z) are contained in the vector
		a=[1, a(1), ... a(p)]
	The input p defines the number of poles in the model.
	The modeling error is returned in err.
	The numerator b(0) is typically set equal to the square
	root of err.
    """
    x = x.flatten()
    N = len(x)
    if p > N:
        print('ERROR: model order too large')
    else:
        X = convm(x, p+1)
        Xq  = X[0:N+p-1,0:p]
        xq1 = -X[1:N+p, 0]
        a = np.linalg.lstsq(Xq, xq1)[0]
        a = np.insert(a, 0, 1)
        err = np.dot(X[0:N+p,0].conj().T, X)
        err = np.dot(err, a)
        err = np.abs(err)
    return a, err


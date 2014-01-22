#!/usr/bin/env python
# vim: set fileencoding=utf-8
""" This file is a Python translation of the MATLAB file prony.m
 
 Python version by RDL 12 Jan 2012
 Copyright notice from prony.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""
 
from __future__ import print_function,division
import sys
import numpy as np
 
from convm import convm
 
def prony(x, p, q):
    """Model a signal using Prony's method
 
    Usage: [b,a,err] = prony(x,p,q)
 
    The input sequence x is modeled as the unit sample response of
    a filter having a system function of the form
        H(z) = B(z)/A(z) 
    The polynomials B(z) and A(z) are formed from the vectors
        b=[b(0), b(1), ... b(q)]
        a=[1   , a(1), ... a(p)]
    The input q defines the number of zeros in the model
    and p defines the number of poles. The modeling error is 
    returned in err.
 
    This comes from Hayes, p. 149, 153, etc
 
    """
    x = x[:]
    N = len(x)
    if p+q >= len(x):
        print('ERROR: model order too large')
        print ("p q len(x) " + str(p) + " " + str(q) + " " + str(len(x)))
        sys.exit(1)
 
    # This formulation uses eq. 4.50, p. 153
    # Set up the convolution matrices
    X = convm(x, p+1)
    Xq = X[q:N+p-1, 0:p]
    xq1 = -X[q+1:N+p,0]
 
    # Solve for denominator coefficients
    if p>0:
        a = np.linalg.lstsq(Xq, xq1)[0]
        a = np.insert(a, 0, 1) # a(0) is 1
    else:
        # all-zero model
        a = np.array(1)
 
    # Solve for the model error
    err = np.dot(x[q+1:N].conj().T,X[q+1:N, 0:p+1])
    err = np.dot(err, a)
 
    # Solve for numerator coefficients
    if q>0:
        # (This is the same as for Pad?)
        b = np.dot(X[0:q+1,0:p+1], a)
    else:
        # all-pole model
        # b(0) is x(0), but a better solution is to match energy
        b = np.sqrt(err)
 
    return (b,a)
 
#function [a,b,err] = prony(x,p,q)
#x   = x(:);
#N   = length(x);
#if p+q>=length(x), error('Model order too large'), end
#X   = convm(x,p+1);
#Xq  = X(q+1:N+p-1,1:p);
#a   = [1;-Xq\X(q+2:N+p,1)];
#b   = X(1:q+1,1:p+1)*a;
#err = x(q+2:N)'*X(q+2:N,1:p+1)*a;
 
def main():
    """Test driver"""
    # From pp. 149-150
    x = np.ones(21)
    p = q = 1
    print('x: {}\np: {}\nq: {}'.format(x,p,q))
    b,a,err = prony(x, p, q)
    print('a: {}\nb: {}\nerr: {}'.format(a,b,err))
 
    # From pp. 152-153
    # Note that these results don't match the book, but they do match the
    # MATLAB version. So I'm either setting things up wrong or this is an
    # errata in the book.
    p = q = 5
    nd = 5
    n = np.arange(11)
    i = np.sinc((n-nd)/2)/2
    b,a,err = prony(i, p, q)
    print('a: {}\nb: {}\nerr: {}'.format(a,b,err))
 
if __name__ == '__main__':
    main()
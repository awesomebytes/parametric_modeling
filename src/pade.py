#!/usr/bin/env python
""" This file is a Python translation of the MATLAB file pade.m

 Python version by RDL 10 Jan 2012
 Copyright notice from pade.m:
 copyright 1996, by M.H. Hayes.  For use with the book 
 "Statistical Digital Signal Processing and Modeling"
 (John Wiley & Sons, 1996).
"""

import sys
import numpy as np

from convm import convm

def pade(x, p, q):
    """Model a signal using the Padé approximation method
    
    Usage: [a,b] = pade(x,p,q)
    
    The input sequence x is modeled as the unit sample response of
    a filter having a system function of the form
        H(z) = B(z)/A(z) 
    The polynomials B(z) and A(z) are formed from the vectors
		b=[b(0), b(1), ... b(q)]
		a=[1   , a(1), ... a(p)]
    The input q defines the number of zeros in the model
    and p defines the number of poles.
    
    This comes from Hayes, p. 138

    """
    if p+q >= len(x):
        print('ERROR: model order too large')
        sys.exit(1)

    # Set up the convolution matrices
    X = convm(x[:], p+1)
    Xq = X[q+1:q+p+1, 1:p+1]
    
    # Solve for the denominator coefficients
    if p>0:    
        a = -np.linalg.lstsq(Xq, x[q+1:q+p+1])[0]
        a = np.insert(a, 0, 1)
    else:
        a = np.array(1)
    
    # Solve for the numerator coefficients
    b = np.dot(X[0:q+1,0:p+1], a)
    
    return (a.flatten(),b.flatten())
    
#function [a,b] = pade(x,p,q)
#x   = x(:);
#X   = convm(x,p+1);
#Xq  = X(q+2:q+p+1,2:p+1);
#a   = [1;-Xq\X(q+2:q+p+1,1)];
#b   = X(1:q+1,1:p+1)*a;


def main():
    """Just a test driver, compare with Hayes pp. 138-140"""
    x = np.array([1, 1.5, 0.75, 0.375, 0.1875, 0.0938])
    pq_array = [(2,0), (0,2), (1,1)]
    for i in xrange(len(pq_array)):
        p,q = pq_array[i]        
        print('For p={} q={}:'.format(p,q))        
        a,b = pade(x, p, q)
        print('a: {}\nb: {}'.format(a,b))
    
if __name__ == '__main__':
    main()


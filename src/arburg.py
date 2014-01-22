#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 20:38 2014

@author: Sammy Pfeiffer
This file pretends to imitate the behaviour of the MATLAB function arburg

This function can't be tested as the main MATLAB functions isn't even working,
gives:

Undefined function 'chkinputdatatype' for input arguments of type 'double'.

But there is still a python version from which we can get it's functionality
http://thomas-cokelaer.info/software/spectrum/html/user/ref_param.html#spectrum.burg.arburg

"""
#import numpy as np
#import scipy
#import matcompat
#from matcompat import *
import spectrum

#def arburg(x, p): # original def
def arburg(args):
    return spectrum.arburg(args[0], args[1], args[2])






##########
#ORIGINAL IMPLEMENTATION AFTER THIS LINES
##########

    # Local Variables: a, E, efp, k, ef, m, varargout, N, p, num, ebp, den, x, eb
    # Function calls: arburg, validateattributes, nargchk, flipud, nargout, issparse, nargin, length, zeros, numel, error, message, conj
    #%ARBURG   AR parameter estimation via Burg method.
    #%   A = ARBURG(X,ORDER) returns the polynomial A corresponding to the AR
    #%   parametric signal model estimate of vector X using Burg's method.
    #%   ORDER is the model order of the AR system.
    #%
    #%   [A,E] = ARBURG(...) returns the final prediction error E (the variance
    #%   estimate of the white noise input to the AR model).
    #%
    #%   [A,E,K] = ARBURG(...) returns the vector K of reflection 
    #%   coefficients (parcor coefficients).
    #%
    #%   % Example:
    #%   %   Estimate input noise variance for AR(4) model.
    #%
    #%   A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
    #%   % Generate noise standard deviations
    #%   % Seed random number generator for reproducible results
    #%   rng default;
    #%   noise_stdz=rand(50,1)+0.5;
    #%   for j=1:50
    #%       y=filter(1,A,noise_stdz(j)*randn(1024,1));
    #%       [ar_coeffs,NoiseVariance(j)]=arburg(y,4);
    #%   end
    #%   %Compare actual vs. estimated variances
    #%   plot(noise_stdz.^2,NoiseVariance,'k*');
    #%   xlabel('Input Noise Variance');
    #%   ylabel('Estimated Noise Variance');
    #%
    #%   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.
    #%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
    #%              Prentice-Hall, 1988, Chapter 7
    #%        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
    #%              Macmillan, 1988, Chapter 5
    #%   Author(s): D. Orofino and R. Losada
    #%   Copyright 1988-2009 The MathWorks, Inc.
    #%   $Revision: 1.12.4.7 $  $Date: 2012/10/29 19:30:37 $
    #matcompat.error(nargchk(2., 2., nargin, 'struct'))
    #% Check the input data type. Single precision is not supported.
    #%try
    #%    chkinputdatatype(x,p);
    #%catch ME
    #%    throwAsCaller(ME);
    #%end
    #validateattributes(x, cellarray(np.hstack(('numeric'))), cellarray(np.hstack(('nonempty', 'finite', 'vector'))), 'arburg', 'X')
    #validateattributes(p, cellarray(np.hstack(('numeric'))), cellarray(np.hstack(('positive', 'integer', 'scalar'))), 'arburg', 'ORDER')
    #if issparse(x):
    #    matcompat.error(message('signal:arburg:Sparse'))
    
    
    #if numel(x)<p+1.:
    #    matcompat.error(message('signal:arburg:InvalidDimension', (p+1.)))
#     x = args[0]
#     p = args[1]
#     
#     x = x.flatten(1)
#     N = length(x)
#     #% Initialization
#     ef = x
#     eb = x
#     a = 1.
#     #% Initial error
#     E = np.dot(x.conj().T, x)/N
#     #% Preallocate 'k' for speed.
#     k = np.zeros(p)
#     for m in np.arange(1., (p)+1):
#         #% Calculate the next order reflection (parcor) coefficient
#         efp = ef[1:]
#         ebp = eb[0:-1]
#         num = np.dot(np.dot(-2., ebp.conj().T), efp)
#         den = np.dot(efp.conj().T, efp)+np.dot(ebp.conj().T, ebp)
#         k[int(m)-1] = num/den
#         #% Update the forward and backward prediction errors
#         ef = efp+np.dot(k[int(m)-1], ebp)
#         eb = ebp+np.dot(k[int(m)-1].conj().T, efp)
#         #% Update the AR coeff.
#         a = np.array(np.vstack((np.hstack((a)), np.hstack((0.)))))+np.dot(k[int(m)-1], np.array(np.vstack((np.hstack((0.)), np.hstack((np.conj(np.flipud(a))))))))
#         #% Update the prediction error
#         E[int((m+1.))-1] = np.dot(1.-np.dot(k[int(m)-1].conj().T, k[int(m)-1]), E[int(m)-1])
#     a = a.flatten(0)
#     #% By convention all polynomials are row vectors
#     varargout = []
#     varargout.append( a)
#     if len(args) >= 2:
#         varargout.append( E[int(0)-1])
#     
#     
#     if len(args) >= 3:
#         varargout.append( k.flatten(1))
#     
#     
#     return [varargout]










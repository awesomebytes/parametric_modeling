#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 20:38 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com
This file pretends to imitate the behaviour of the MATLAB function aryule

Using spectrum implementation:
http://thomas-cokelaer.info/software/spectrum/html/user/ref_param.html#spectrum.yulewalker.aryule

"""
import numpy as np
import spectrum

def aryule(x, p):
    [A, E, K] = spectrum.aryule(x, p)
    A = np.hstack((1, A)) # MATLAB adds the first "1.0"
    return A, E, K

    # Local Variables: a, e, k, nx, p, R, x, mx
    # Function calls: aryule, nargchk, min, issparse, nargin, length, isempty, error, levinson, message, xcorr, round, size
    #%ARYULE   AR parameter estimation via Yule-Walker method.
    #%   A = ARYULE(X,ORDER) returns the polynomial A corresponding to the AR
    #%   parametric signal model estimate of vector X using the Yule-Walker
    #%   (autocorrelation) method.  ORDER is the model order of the AR system. 
    #%   This method solves the Yule-Walker equations by means of the Levinson-
    #%   Durbin recursion.
    #%
    #%   [A,E] = ARYULE(...) returns the final prediction error E (the variance
    #%   estimate of the white noise input to the AR model).
    #%
    #%   [A,E,K] = ARYULE(...) returns the vector K of reflection coefficients.
    #%   
    #%   % Example:
    #%   %   Estimate model order using decay of reflection coefficients.
    #%
    #%   rng default;
    #%   y=filter(1,[1 -0.75 0.5],0.2*randn(1024,1));
    #%
    #%   % Create AR(2) process
    #%   [ar_coeffs,NoiseVariance,reflect_coeffs]=aryule(y,10);
    #%
    #%   % Fit AR(10) model
    #%   stem(reflect_coeffs); axis([-0.05 10.5 -1 1]);
    #%   title('Reflection Coefficients by Lag'); xlabel('Lag');
    #%   ylabel('Reflection Coefficent');
    #%
    #%   See also PYULEAR, ARMCOV, ARBURG, ARCOV, LPC, PRONY.
    #%   Ref: S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
    #%              Macmillan, 1988, Chapter 5
    #%        M. Hayes, STATISTICAL DIGITAL SIGNAL PROCESSING AND MODELING, 
    #%              John Wiley & Sons, 1996, Chapter 8
    #%   Author(s): R. Losada
    #%   Copyright 1988-2004 The MathWorks, Inc.
    #%   $Revision: 1.12.4.6 $  $Date: 2012/10/29 19:30:38 $
#     matcompat.error(nargchk(2., 2., nargin, 'struct'))
#     #% Check the input data type. Single precision is not supported.
#     #%try
#     #%    chkinputdatatype(x,p);
#     #%catch ME
#     #%    throwAsCaller(ME);
#     #%end
#     [mx, nx] = matcompat.size(x)
#     if isempty(x) or length(x)<p or matcompat.max(mx, nx) > 1.:
#         matcompat.error(message('signal:aryule:InvalidDimensions'))
#     elif isempty(p) or not p == np.round(p):
#         matcompat.error(message('signal:aryule:MustBeInteger'))
#         
#     
#     if issparse(x):
#         matcompat.error(message('signal:aryule:Sparse'))
#     
#     
#     R = plt.xcorr(x, p, 'biased')
#     [a, e, k] = levinson(R[int(p+1.)-1:], p)
#     return [a, e, k]
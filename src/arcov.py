#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 20:38 2014

@author: Sammy Pfeiffer
This file pretends to imitate the behaviour of the MATLAB function arcov

Using spectrum implementation:
http://thomas-cokelaer.info/software/spectrum/html/user/ref_psd_other.html#spectrum.covar.arcovar
"""
import numpy as np
import spectrum

def arcov(x, p):
    [A, E] = spectrum.covar.arcovar(x,p)
    A = np.hstack((1,A)) # MATLAB gives back initial value 1, so we do the same 
    return A, E

    # Local Variables: a, msgobj, p, msg, x, e
    # Function calls: nargchk, arcov, nargin, isempty, error, arparest
    #%ARCOV   AR parameter estimation via covariance method.
    #%   A = ARCOV(X,ORDER) returns the polynomial A corresponding to the AR
    #%   parametric signal model estimate of vector X using the Covariance method.
    #%   ORDER is the model order of the AR system.
    #%
    #%   [A,E] = ARCOV(...) returns the variance estimate E of the white noise
    #%   input to the AR model.
    #%
    #%   See also PCOV, ARMCOV, ARBURG, ARYULE, LPC, PRONY.
    #%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
    #%              Prentice-Hall, 1988, Chapter 7
    #%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
    #%              Prentice-Hall, 1997, Chapter 3
    #%   Author(s): R. Losada and P. Pacheco
    #%   Copyright 1988-2002 The MathWorks, Inc.
    #%   $Revision: 1.13.4.3 $  $Date: 2011/05/13 18:06:51 $
    #matcompat.error(nargchk(2., 2., nargin, 'struct'))
#     [a, e, msg, msgobj] = arparest(x, p, 'covariance') 
#     if not isempty(msg):
#         matcompat.error(msgobj)
#     
#     
#     #% [EOF] - arcov.m
#     return [a, e]
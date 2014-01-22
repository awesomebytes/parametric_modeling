#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jan 22 20:38 2014

@author: Sammy Pfeiffer
This file pretends to imitate the behaviour of the MATLAB function armcov

Using spectrum implementation:
http://thomas-cokelaer.info/software/spectrum/html/user/ref_psd_other.html
search for modcovar (can't direct link it)
"""
import numpy as np
# import scipy
# from matcompat import *
import spectrum

def armcov(x, p):
    [E, A] = spectrum.modcovar(x, int(p)) # We compute this one because gives back same number of elements in E
    number_of_elements = len(E)
    [E, A, ISTAT] =spectrum.modcovar_marple(x, int(p)) # works slower but is more accurate with the error than modcovar
    E = E[:number_of_elements]
    E = np.hstack((1, E))
    
    return E, A

### ORIGINAL IMPLEMENTATION
    # Local Variables: a, msgobj, p, msg, x, e
    # Function calls: nargchk, nargin, isempty, error, arparest, armcov
    #%ARMCOV   AR parameter estimation via modified covariance method.
    #%   A = ARMCOV(X,ORDER) returns the polynomial A corresponding to the AR
    #%   parametric signal model estimate of vector X using the Modified Covariance
    #%   method. ORDER is the model order of the AR system. 
    #%
    #%   [A,E] = ARMCOV(...) returns the variance estimate E of the white noise
    #%   input to the AR model.
    #%
    #%   See also PMCOV, ARCOV, ARBURG, ARYULE, LPC, PRONY.
    #%   References:
    #%     [1] S. Lawrence Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATIONS,
    #%              Prentice-Hall, 1987, Chapter 8
    #%     [2] Steven M. Kay, MODERN SPECTRAL ESTIMATION THEORY & APPLICATION,
    #%              Prentice-Hall, 1988, Chapter 7
    #%   Author(s): R. Losada and P. Pacheco
    #%   Copyright 1988-2002 The MathWorks, Inc.
    #%   $Revision: 1.13.4.3 $  $Date: 2011/05/13 18:06:53 $
#     matcompat.error(nargchk(2., 2., nargin, 'struct'))
#     [a, e, msg, msgobj] = arparest(x, p, 'modified')
#     if not isempty(msg):
#         matcompat.error(msgobj)
#     
#     
#     #% [EOF] - armcov.m
#     return [a, e]
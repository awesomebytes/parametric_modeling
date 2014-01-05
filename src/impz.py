#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 5 17:48:01 2014

@author: Sammy Pfeiffer
This file contains a python implementation of impz matlab function
"""
from scipy import zeros
from scipy.signal import lfilter

def impz(b,a):
    """Pseudo implementation of the impz method of MATLAB"""
#% Compute time vector
# M = 0;  NN = [];
# if isempty(N)
#   % if not specified, determine the length
#   if isTF
#     N = impzlength(b,a,.00005);
#   else
#     N  = impzlength(b,.00005);
#   end
    N = len(b) * len(b) * len(b) # MATLAB AUTOFINDS THE SIZE HERE... 
    x = zeros(N)
    x[0] = 1
    h = lfilter(b,a, x)
    return h
    
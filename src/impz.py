#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 5 17:48:01 2014

@author: Sammy Pfeiffer
This file contains a python implementation of impz matlab function
"""
from scipy import zeros
from scipy.signal import lfilter
import numpy as np

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
    p = np.roots(a)
    N = stableNmarginal_length(p, 0.00005, 0)
    N = len(b) * len(b) * len(b) # MATLAB AUTOFINDS THE SIZE HERE... 
    #TODO: Implement some way of finding the autosieze of this... I used a couple of examples... matlab gave 43 as length we give 64
    x = zeros(N)
    x[0] = 1
    h = lfilter(b,a, x)
    return h
    
    
def stableNmarginal_length(p, tol, delay):
    """ % Determine the length for an unstable filter
    %minimum height is .00005 original amplitude:"""
    
# ind = find(abs(p-1)<1e-5); # does nothing in our example case
# p(ind) = -p(ind);    % treat constant as Nyquist # does nothing
# ind = find(abs(abs(p)-1)<1e-5);     # does nothing too...  
# periods = 5*max(2*pi./abs(angle(p(ind)))); % five periods
# p(ind) = [];   % get rid of unit circle poles
# [maxp,maxind] = max(abs(p));
# if isempty(p)   % pure oscillator
#     N = periods;
# elseif isempty(ind)   % no oscillation
#     N = mltplcty(p,maxind)*log10(tol)/log10(maxp) + delay;
# else    % some of both
#     N = max(periods, ...
#         mltplcty(p,maxind)*log10(tol)/log10(maxp) ) + delay;
    return 1

# function m = mltplcty( p, ind, tol)
# %MLTPLCTY  Multiplicity of a pole
# %   MLTPLCTY(P,IND,TOL) finds the multiplicity of P(IND) in the vector P
# %   with a tolerance of TOL.  TOL defaults to .001.
# 
# if nargin<3
#     tol = .001;
# end
# 
# [mults,indx]=mpoles(p,tol);
# 
# m = mults(indx(ind));
# for i=indx(ind)+1:length(mults)
#     if mults(i)>m
#         m = m + 1;
#     else
#         break;
#     end
# end
    
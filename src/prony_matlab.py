#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 6 20:54:46 2014

@author: Sammy Pfeiffer
This file pretends to imitate the behaviour of the MATLAB function prony
"""

from scipy.linalg import toeplitz
import numpy as np

# function [b,a] = prony(h, nb ,na)
def prony_matlab(h, nb, na):
# %PRONY Prony's method for time-domain IIR filter design.
# %   [B,A] = PRONY(H, NB, NA) finds a filter with numerator order
# %   NB, denominator order NA, and having the impulse response in
# %   vector H.   The IIR filter coefficients are returned in
# %   length NB+1 and NA+1 row vectors B and A, ordered in
# %   descending powers of Z.  H may be real or complex.
# %
# %   If the largest order specified is greater than the length of H,
# %   H is padded with zeros.
# %
# %   % Example:
# %   %   Fit an IIR model to an impulse response of a lowpass filter.
# %
# %   [b,a] = butter(4,0.2);
# %   impulseResp = impz(b,a);                % obtain impulse response
# %   denOrder=4; numOrder=4;                 % system function of order 4
# %   [Num,Den]=prony(impulseResp,numOrder,denOrder);
# %   subplot(211);                           % impulse response and input
# %   stem(impz(Num,Den,length(impulseResp)));   
# %   title('Impulse Response with Prony Design');
# %   subplot(212);
# %   stem(impulseResp); title('Input Impulse Response');
# %
# %   See also STMCB, LPC, BUTTER, CHEBY1, CHEBY2, ELLIP, INVFREQZ.
# 
# %   Author(s): L. Shure, 5-17-88
# %              L. Shure, 12-17-90, revised
# %   Copyright 1988-2012 The MathWorks, Inc.
# %   $Revision: 1.7.4.1.2.1 $  $Date: 2013/01/02 17:47:48 $
# 
# %   References:
# %     [1] T.W. Parks and C.S. Burrus, Digital Filter Design,
# %         John Wiley and Sons, 1987, p226.
# 
# K = length(h) - 1;
# M = nb; N = na;
    K = len(h) - 1
    M = nb
    N = na

# if K <= max(M,N)      % zero-pad input if necessary
#     K = max(M,N)+1;
#     h(K+1) = 0;
# end
    if K <= max(M,N):
        K = max(M,N) + 1
        h[K+1] = 0 # probable problem with indices!
# c = h(1);
    c = h[0] # probable problem with indices!
# if c==0    % avoid divide by zero
#     c=1;
# end
    if c == 0:
        c = 1  
# H = toeplitz(h/c,[1 zeros(1,K)]);
    print "toeplitz second part:"
    print [1, np.zeros(1,K)]
    H = toeplitz(h/c, [1, np.zeros(1,K)]) # probable problem with indices!
    print "H (thanks to toeplitz) is:"
    print H
# % K+1 by N+1
# if (K > N)
#     H(:,(N+2):(K+1)) = [];
# end
    if K > N:
        H[:,N+2:K+1] = [] # probable problem with indices!

# % Partition H matrix
# H1 = H(1:(M+1),:);    % M+1 by N+1
    H1 = H[0:M+1,:] # probable problem with indices!
# h1 = H((M+2):(K+1),1);    % K-M by 1
    h1 = H[M+2:K+1,0] # probable problem with indices!
# H2 = H((M+2):(K+1),2:(N+1));    % K-M by N
    H2 = H[M+2:K+1,1:N+1] # probable problem with indices!
# a = [1; -H2\h1].';
    #a = [0, -H2\h1]
    #\
    # Matrix left division
    # x = A\B is the solution to the equation Ax = B. Matrices A and B must have the same number of rows.
# b = c*a*H1.';
    #b = c*a*H1
    b = 0
    a = 0
    return [b, a]
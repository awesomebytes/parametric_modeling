#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 6 20:54:46 2014
Redone on Jan 22 12:00:00 2014

@author: Sammy Pfeiffer
@email: sammypfeiffer@gmail.com
This file pretends to imitate the behavior of the MATLAB function prony
"""
import numpy as np
from scipy.linalg import toeplitz


def prony(h, nb, na):
    """From MATLAB:
    %PRONY Prony's method for time-domain IIR filter design.
    %   [B,A] = PRONY(H, NB, NA) finds a filter with numerator order
    %   NB, denominator order NA, and having the impulse response in
    %   vector H.   The IIR filter coefficients are returned in
    %   length NB+1 and NA+1 row vectors B and A, ordered in
    %   descending powers of Z.  H may be real or complex.
    %
    %   If the largest order specified is greater than the length of H,
    %   H is padded with zeros.
    """

    # Local Variables: a, c, b, h1, h, nb, H1, M, N, H2, H, na, H2_minus, K
    # Function calls: max, length, prony, zeros, toeplitz
    #%PRONY Prony's method for time-domain IIR filter design.
    #%   [B,A] = PRONY(H, NB, NA) finds a filter with numerator order
    #%   NB, denominator order NA, and having the impulse response in
    #%   vector H.   The IIR filter coefficients are returned in
    #%   length NB+1 and NA+1 row vectors B and A, ordered in
    #%   descending powers of Z.  H may be real or complex.
    #%
    #%   If the largest order specified is greater than the length of H,
    #%   H is padded with zeros.
    #%
    #%   % Example:
    #%   %   Fit an IIR model to an impulse response of a lowpass filter.
    #%
    #%   [b,a] = butter(4,0.2);
    #%   impulseResp = impz(b,a);                % obtain impulse response
    #%   denOrder=4; numOrder=4;                 % system function of order 4
    #%   [Num,Den]=prony(impulseResp,numOrder,denOrder);
    #%   subplot(211);                           % impulse response and input
    #%   stem(impz(Num,Den,length(impulseResp)));   
    #%   title('Impulse Response with Prony Design');
    #%   subplot(212);
    #%   stem(impulseResp); title('Input Impulse Response');
    #%
    #%   See also STMCB, LPC, BUTTER, CHEBY1, CHEBY2, ELLIP, INVFREQZ.
    #%   Author(s): L. Shure, 47-88
    #%              L. Shure, 117-90, revised
    #%   Copyright 1988-2012 The MathWorks, Inc.
    #%   $Revision: 1.7.4.1.2.1 $  $Date: 2013/01/02 17:47:48 $
    #%   References:
    #%     [1] T.W. Parks and C.S. Burrus, Digital Filter Design,
    #%         John Wiley and Sons, 1987, p226.
    K = len(h)-1
    M = nb
    N = na
    if K<=max(M, N):
        #% zero-pad input if necessary
        K = max(M, N)+1
        h[K+1] = 0
    
    c = h[0]
    if c == 0:
        #% avoid divide by zero
        c = 1
    
#     second_part_toeplitz = np.zeros(K+1)
#     second_part_toeplitz[0] = 1
    H = toeplitz(h / c, np.array(np.hstack((1, np.zeros(K)))))
    #% K+1 by N+1
    if K > N:
        #H[:,int(N+2)-1:K+1] = np.array([])
        H = H[:,0:N+1] # Here we are just getting rid of all the columns after N+1
    
    
    #% Partition H matrix
    H1 = H[0:M+1,:]
    #% M+1 by N+1
    h1 = H[M+1:K+1,0]
    #% K-M by 1
    H2 = H[M:K,0:N]
    #% K-M by N
    H2_minus = -H2
    a_right = np.linalg.lstsq(-H2,h1)[0]
    a_left = np.array([1])
    a = np.append(a_left, a_right)
    #a = np.array(np.vstack((np.hstack((1)), np.hstack((solve(H2_minus, h1)))))).T
    b = np.dot(np.dot(c, a), H1.T)
    return [b, a]
#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:30:46 2013

@author: Sammy Pfeiffer
This file pretends to imitate the behaviour of the MATLAB function with the same name.
"""

from prony_matlab import prony_matlab
from convmtx import convmtx
from scipy.signal import lfilter
import numpy as np

def stmcb(x, u_in=None, q=None, p=None, niter=5, a_in=None):
    """
    %STMCB Compute linear model via Steiglitz-McBride iteration
   [B,A] = stmcb(H,NB,NA) finds the coefficients of the system 
   B(z)/A(z) with approximate impulse response H, NA poles and 
   NB zeros.

   [B,A] = stmcb(H,NB,NA,N) uses N iterations.  N defaults to 5.

   [B,A] = stmcb(H,NB,NA,N,Ai) uses the vector Ai as the initial 
   guess at the denominator coefficients.  If you don't specify Ai, 
   STMCB uses [B,Ai] = PRONY(H,0,NA) as the initial conditions.

   [B,A] = STMCB(Y,X,NB,NA,N,Ai) finds the system coefficients B and 
   A of the system which, given X as input, has Y as output.  N and Ai
   are again optional with default values of N = 5, [B,Ai] = PRONY(Y,0,NA).
   Y and X must be the same length.

   % Example:
   %   Approximate the impulse response of a Butterworth filter with a 
   %   system of lower order.

   [b,a] = butter(6,0.2);              % Butterworth filter design
   h = filter(b,a,[1 zeros(1,100)]);   % Filter data using above filter
   freqz(b,a,128)                      % Frequency response 
   [bb,aa] = stmcb(h,4,4);                  
   figure; freqz(bb,aa,128)

   See also PRONY, LEVINSON, LPC, ARYULE.

   Author(s): Jim McClellan, 2-89
              T. Krauss, 4-22-93, new help and options
   Copyright 1988-2004 The MathWorks, Inc.
   $Revision: 1.8.4.6 $  $Date: 2012/10/29 19:32:10 $
    """
    
    #error(nargchk(3,6,nargin,'struct'))
    
#TODO: fit the full definition of the function
#   if length(u_in) == 1,
#     if nargin == 3,
#         niter = 5; p = q; q = u_in;
#         a_in = prony(x,0,p);
#     elseif nargin == 4,
#         niter = p; p = q; q = u_in; 
#         a_in = prony(x,0,p);
#     elseif nargin == 5,
#         a_in = niter; niter = p; p = q; q = u_in; 
#     end
#     u_in = zeros(size(x));
#     u_in(1) = 1;         % make a unit impulse whose length is same as x
# else
#     if length(u_in)~=length(x),
#         error(message('signal:stmcb:InvalidDimensions'))
#     end
    if len(u_in) != len(x):
        print "stmbc:"
        print "Invalid dimensions on u_in and x, must be of the same size: " + str(len(u_in)) + "!=" + str(len(x)) 
        exit(0)
#     if nargin < 6
#        [b,a_in] = prony(x,0,p);
#     end
    if a_in == None:
        [b, a_in] = prony_matlab(x, 0, p)
        print "p is:"
        print p
        print "b is: "
        print b
        print "a_in is:"
        print a_in
#     if nargin < 5
#        niter = 5;
#     end
    # nargin already initialized as 5 on function definition, check not necessary
# end
    
    
    
# a = a_in;
# N = length(x);
    a = a_in
    N = len(x)
# for i=1:niter
    for i in range(niter):
#    u = filter( 1, a, x );
        print "a is: " + str(a)
        ## MATLAB
        # a = 
        #  1.0000
        # -1.9803
        #  0.9806
        # python a is: 1

        #print "x is: " + str(x)
        u = lfilter(1, a, x)
        print u
        exit(0)
#    v = filter( 1, a, u_in );
        v = lfilter(1, a, u_in)
#    C1 = convmtx(u(:),p+1);
        C1 = convmtx(u,p+1)
#    C2 = convmtx(v(:),q+1);
        C2 = convmtx(v,q+1)
#    T = [ -C1(1:N,:) C2(1:N,:) ];
        T = np.array([-C1[0:N,:],
                       C2[0:N,:]])
#    c = T(:,2:p+q+2)\(-T(:,1));   % move 1st column to RHS and do least-squares
# Si la matriz no es cuadrada: numpy.linalg.lstsq
# Si la matriz es cuadrada: numpy.linalg.solve
        if T.shape[0] != T.shape[1]:
            c = np.linalg.lstsq(a, b)

#    a = [1; c(1:p)];                % denominator coefficients
#    b = c(p+1:p+q+1);               % numerator coefficients
# end
# a=a.';
# b=b.';
    a = b = 0
    return b, a
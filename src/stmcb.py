
import numpy as np
import scipy
from scipy.linalg import solve, lstsq
from scipy.signal import lfilter
import matcompat
from matcompat import *
from convmtx import convmtx
from prony import prony


# def printInfo(name, thing_to_print):
#     print name + " is: ", thing_to_print
#     print name + " shape is:", thing_to_print.shape

#def stmcb(x, u_in=None, q=None, p=None, niter=5, a_in=None): # Old definition
#function [b,a] = stmcb_test( x, u_in, q, p, niter, a_in ) # matlab definition
def stmcb(*args):

    # Local Variables: T_sub2, T_sub1, T_minus, C2, C1, a_in, N, u_in, T, a, niter, c, b, i, q, p, u, v, x, T_left, T_right, C1_minus
    # Function calls: convmtx, filter, prony, nargchk, stmcb, nargin, length, zeros, error, message, size
    #%STMCB Compute linear model via Steiglitz-McBride iteration
    #%   [B,A] = stmcb(H,NB,NA) finds the coefficients of the system 
    #%   B(z)/A(z) with approximate impulse response H, NA poles and 
    #%   NB zeros.
    #%
    #%   [B,A] = stmcb(H,NB,NA,N) uses N iterations.  N defaults to 5.
    #%
    #%   [B,A] = stmcb(H,NB,NA,N,Ai) uses the vector Ai as the initial 
    #%   guess at the denominator coefficients.  If you don't specify Ai, 
    #%   STMCB uses [B,Ai] = PRONY(H,0,NA) as the initial conditions.
    #%
    #%   [B,A] = STMCB(Y,X,NB,NA,N,Ai) finds the system coefficients B and 
    #%   A of the system which, given X as input, has Y as output.  N and Ai
    #%   are again optional with default values of N = 5, [B,Ai] = PRONY(Y,0,NA).
    #%   Y and X must be the same length.
    #%
    #%   % Example:
    #%   %   Approximate the impulse response of a Butterworth filter with a 
    #%   %   system of lower order.
    #%
    #%   [b,a] = butter(6,0.2);              % Butterworth filter design
    #%   h = filter(b,a,[1 zeros(1,100)]);   % Filter data using above filter
    #%   freqz(b,a,128)                      % Frequency response 
    #%   [bb,aa] = stmcb(h,4,4);                  
    #%   figure; freqz(bb,aa,128)
    #%
    #%   See also PRONY, LEVINSON, LPC, ARYULE.
    #%   Author(s): Jim McClellan, 2-89
    #%   	       T. Krauss, 4-22-93, new help and options
    #%   Copyright 1988-2004 The MathWorks, Inc.
    #%   $Revision: 1.8.4.6 $  $Date: 2012/10/29 19:32:10 $
    #matcompat.error(nargchk(3, 6, nargin, 'struct'))
    x = args[0]
    # Default assignments
    u_in = args[1]
    q = args[2]
    p = args[3]
    if type(u_in) == type(1):
        if len(args) == 3:
            niter = 5
            p = args[2]
            q = args[1]
            [a_in, tmp] = prony(x, 0, p) # In Python there is no way to just get one of the two items returned
            
        elif len(args) == 4:
            niter = args[3]
            p = args[2]
            q = u_in
            [a_in, tmp] = prony(x, 0, p)
            
        elif len(args) == 5:
            a_in = args[4]
            niter = args[3]
            p = args[2]
            q = u_in
            
        u_in = np.zeros(matcompat.size(x))
        u_in[0] = 1.
        #% make a unit impulse whose length is same as x
    else:
        if len(u_in) != len(x):
            print "stmcb:"
            print "Invalid dimensions on u_in and x, must be of the same size: " + str(len(u_in)) + "!=" + str(len(x)) 
            exit(0)
            
        if len(args) < 6:
            [b, a_in] = prony(x, 0, args[3])
        
        if len(args) < 5:
            niter = 5
        

    a = a_in
    N = len(x)
    for i in range(niter):
        u = lfilter([1], a, x)
        v = lfilter([1], a, u_in)
        C1 = convmtx(u, (p+1))
        C2 = convmtx(v, (q+1))
        T_left = -C1[0:N,:]
        T_right = C2[0:N,:]
        T = np.hstack((T_left, T_right))
        #%T = [ C1_minus(1:N,:) C2(1:N,:) ];
        T_minus = -T
        T_left = T[:,1:p+q+2]
        T_right = T_minus[:,0]
        # If not squared matrix: numpy.linalg.lstsq
        # If squared matrix: numpy.linalg.solve
        if T.shape[0] != T.shape[1]:
            [c, residuals, rank, singular_values] = lstsq(T_left, T_right) # lstsq in python returns more stuff
        else:
            c = solve(T_left, T_right)
        #%    c = T(:,2:p+q+2)\(T_minus(:,1));   % move 1st column to RHS and do least-squares
        a_left = np.array([1])
        a_right = c[:p]
        a = np.hstack((a_left, a_right))
        #% denominator coefficients
        b = c[p:p+q+1]
        #% numerator coefficients
    a = a.T
    b = b.T
    return [b, a]
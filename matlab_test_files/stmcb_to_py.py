
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def stmcb(x, u_in, q, p, niter, a_in):

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
    matcompat.error(nargchk(3., 6., nargin, 'struct'))
    if length(u_in) == 1.:
        if nargin == 3.:
            niter = 5.
            p = q
            q = u_in
            a_in = prony(x, 0., p)
        elif nargin == 4.:
            niter = p
            p = q
            q = u_in
            a_in = prony(x, 0., p)
            
        elif nargin == 5.:
            a_in = niter
            niter = p
            p = q
            q = u_in
            
        
        u_in = np.zeros(matcompat.size(x))
        u_in[0] = 1.
        #% make a unit impulse whose length is same as x
    else:
        if length(u_in) != length(x):
            matcompat.error(message('signal:stmcb:InvalidDimensions'))
        
        
        if nargin<6.:
            [b, a_in] = prony(x, 0., p)
        
        
        if nargin<5.:
            niter = 5.
        
        
        
    
    a = a_in
    N = length(x)
    for i in np.arange(1., (niter)+1):
        u = filter(1., a, x)
        v = filter(1., a, u_in)
        C1 = convmtx(u.flatten(1), (p+1.))
        C2 = convmtx(v.flatten(1), (q+1.))
        C1_minus = -C1
        T_sub1 = C1_minus[0:N,:]
        T_sub2 = C2[0:N,:]
        T = np.array(np.hstack((T_sub1, T_sub2)))
        #%T = [ C1_minus(1:N,:) C2(1:N,:) ];
        T_minus = -T
        T_left = T[:,1:p+q+2.]
        T_right = T_minus[:,0]
        c = linalg.solve(T_left, T_right)
        #%    c = T(:,2:p+q+2)\(T_minus(:,1));   % move 1st column to RHS and do least-squares
        a = np.array(np.vstack((np.hstack((1.)), np.hstack((c[0:p])))))
        #% denominator coefficients
        b = c[int(p+1.)-1:p+q+1.]
        #% numerator coefficients
        
    a = a.T
    b = b.T
    return [b, a]
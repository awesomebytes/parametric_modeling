
import numpy as np
import scipy
from matcompat import *

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def invfreqs(g, w, varargin):

    # Local Variables: realFlag, cg, realStr, D31, gndir, t1, cw, nk, kom, nm, na, nb, Vcap, V1, rg, pf, rw, tol, maxiter, varargin, cwf, wf, ll, D, rwf, D32, Dva, Dvb, gaussFlag, verb, GC, e, th, a, OM, b, Vd, g, k, l, st, indg, R, inda, t, w, indb, D3
    # Function calls: disp, all, invfreqs, deal, ischar, int2str, warning, apolystab, home, message, size, getString, sqrt, clc, zeros, norm, real, nargchk, max, nargin, ones, isempty, lower, length, num2str, error, strcmp
    #%INVFREQS Analog filter least squares fit to frequency response data.
    #%   [B,A] = INVFREQS(H,W,nb,na) gives real numerator and denominator 
    #%   coefficients B and A of orders nb and na respectively, where
    #%   H is the desired complex frequency response of the system at frequency
    #%   points W, and W contains the frequency values in radians/s.
    #%   INVFREQS yields a filter with real coefficients.  This means that it is 
    #%   sufficient to specify positive frequencies only; the filter fits the data 
    #%   conj(H) at -W, ensuring the proper frequency domain symmetry for a real 
    #%   filter.
    #%
    #%   [B,A]=INVFREQS(H,W,nb,na,Wt) allows the fit-errors to the weighted
    #%   versus frequency.  LENGTH(Wt)=LENGTH(W)=LENGTH(H).
    #%   Determined by minimization of sum |B-H*A|^2*Wt over the freqs in W.
    #%
    #%   [B,A] = INVFREQS(H,W,nb,na,Wt,ITER) does another type of fit:
    #%   Sum |B/A-H|^2*Wt is minimized with respect to the coefficients in B and
    #%   A by numerical search in at most ITER iterations.  The A-polynomial is 
    #%   then constrained to be stable.  [B,A]=INVFREQS(H,W,nb,na,Wt,ITER,TOL)
    #%   stops the iterations when the norm of the gradient is less than TOL.
    #%   The default value of TOL is 0.01.  The default value of Wt is all ones.
    #%   This default value is also obtained by Wt=[].
    #%
    #%   [B,A]=INVFREQS(H,W,nb,na,Wt,ITER,TOL,'trace') provides a textual
    #%   progress report of the iteration.
    #%
    #%   [B,A] = INVFREQS(H,W,'complex',NB,NA,...) creates a complex filter.  In 
    #%   this case, no symmetry is enforced. 
    #%
    #%   % Example:
    #%   %   Convert a simple transfer function to frequency response data and
    #%   %   then back to the original filter coefficients. If the system is
    #%   %   unstable, use invfreqs's iterative algorithm to find a stable 
    #%   %   approximation to the system.
    #%
    #%   b = [1 2 3 2 3];                % Numerator coefficients
    #%   a = [1 2 3 2 1 4];              % Denominator coefficients
    #%   [h,w] = freqs(b,a,64);
    #%   [bb,aa] = invfreqs(h,w,4,5) %  aa has poles in the right half-plane.
    #%   fprintf('Stable Approximation to the system:')
    #%   [bbb,aaa] = invfreqs(h,w,4,5,[],30) % stable approximation to system
    #%
    #%   See also FREQZ, FREQS, INVFREQZ.
    #%   Author(s): J.O. Smith and J.N. Little, 4-23-86
    #%              J.N. Little, 4-27-88, revised
    #%              Lennart Ljung, 9-21-92, rewritten
    #%              T. Krauss, 10-22-92, trace mode made optional
    #%   Copyright 1988-2011 The MathWorks, Inc.
    #%     
    #% calling sequence is
    #%function [b,a]=invfreqs(g,w,nb,na,wf,maxiter,tol,pf)
    #% OR
    #%function [b,a]=invfreqs(g,w,'complex',nb,na,wf,maxiter,tol,pf)
    matcompat.error(nargchk(4., 9., nargin, 'struct'))
    if ischar(varargin.cell[0]):
        realStr = lower(varargin.cell[0])
        varargin[0] = np.array([])
    else:
        realStr = 'real'
        
    
    gaussFlag = length(varargin) > 3.
    #% run Gauss-Newton algorithm or not?
    if length(varargin)<6.:
        varargin.cell[5] = np.array([])
        #% pad varargin with []'s
    
    
    [nb, na, wf, maxiter, tol, pf] = deal(varargin.cell[:])
    _switch_val=realStr
    if False: # switch 
        pass
    elif _switch_val == 'real':
        realFlag = 1.
    elif _switch_val == 'complex':
        realFlag = 0.
    else:
        matcompat.warning(message('signal:invfreqs:InvalidParam', realStr))
        realFlag = 0.
    
    nk = 0.
    #% The code is prepared for constraining the numerator to 
    #% begin with nk zeros.
    nb = nb+nk+1.
    if isempty(pf):
        verb = 0.
    elif strcmp(pf, 'trace'):
        verb = 1.
        
    else:
        matcompat.error(message('signal:invfreqs:NotSupported', pf))
        
    
    if isempty(wf):
        wf = np.ones(length(w), 1.)
    
    
    wf = np.sqrt(wf)
    if length(g) != length(w):
        matcompat.error(message('signal:invfreqs:UnmatchedLengths', 'H', 'W'))
    
    
    if length(wf) != length(w):
        matcompat.error(message('signal:invfreqs:UnmatchedLengths', 'Wt', 'W'))
    
    
    #% if any( w(:)<0 ) && realFlag
    #%   warning(message('signal:invfreqs:InvalidWParam', 'W', 'INVFREQS', 'complex'))
    #% end
    [rw, cw] = matcompat.size(w)
    if rw > cw:
        w = w.conj().T
    
    
    [rg, cg] = matcompat.size(g)
    if cg > rg:
        g = g.T
    
    
    [rwf, cwf] = matcompat.size(wf)
    if cwf > rwf:
        wf = wf.conj().T
    
    
    nm = matcompat.max((na+1.), (nb+nk))
    indb = np.arange(nb, (1.)+(-1.), -1.)
    indg = np.arange(na+1., (1.)+(-1.), -1.)
    inda = np.arange(na, (1.)+(-1.), -1.)
    OM = np.ones(1., length(w))
    for kom in np.arange(1., (nm-1.)+1):
        OM = np.array(np.vstack((np.hstack((OM)), np.hstack(((1i.*w)**kom)))))
        
    #%
    #% Estimation in the least squares case:
    #%
    Dva = OM[int(inda)-1,:].T*np.dot(g, np.ones(1., na))
    Dvb = -OM[int(indb)-1,:].T
    D = np.array(np.hstack((Dva, Dvb)))*np.dot(wf, np.ones(1., (na+nb)))
    R = np.dot(D.conj().T, D)
    Vd = np.dot(D.conj().T, -g*OM[int((na+1.))-1,:].T*wf)
    if realFlag:
        R = np.real(R)
        Vd = np.real(Vd)
    
    
    th = linalg.solve(R, Vd)
    a = np.array(np.hstack((1., th[0:na].T)))
    b = np.array(np.hstack((np.zeros(1., nk), th[int(na+1.)-1:na+nb].T)))
    if not gaussFlag:
        return []
    
    
    #% Now for the iterative minimization
    if isempty(maxiter):
        maxiter = 30.
    
    
    if isempty(tol):
        tol = 0.01
    
    
    #% Stabilizing the denominator:
    a = apolystab(a, realFlag)
    #% The initial estimate:
    GC = (np.dot(b, OM[int(indb)-1,:])/np.dot(a, OM[int(indg)-1,:])).T
    e = (GC-g)*wf
    Vcap = np.dot(e.conj().T, e)
    t = np.array(np.hstack((a[1:na+1.], b[int(nk+1.)-1:nk+nb]))).T
    if verb:
        #% invfreqz using same messages
    clc
    np.disp(np.array(np.hstack(('  ', getString(message('signal:invfreqs:INITIALESTIMATE'))))))
    np.disp(np.array(np.hstack((getString(message('signal:invfreqs:CurrentFit')), num2str(Vcap)))))
    np.disp(getString(message('signal:invfreqs:Parvector')))
    np.disp(t)
    
    #% %
    #% ** the minimization loop **
    #%
    gndir = 2.*tol+1.
    l = 0.
    st = 0.
    while np.all(np.array(np.hstack((linalg.norm(gndir) > tol, l<maxiter, st != 1.)))):
        l = l+1.
        #%     * compute gradient *
        D31 = OM[int(inda)-1,:].T*np.dot(-GC/np.dot(a, OM[int(indg)-1,:]).T, np.ones(1., na))
        D32 = OM[int(indb)-1,:].T/np.dot(np.dot(a, OM[int(indg)-1,:]).T, np.ones(1., nb))
        D3 = np.array(np.hstack((D31, D32)))*np.dot(wf, np.ones(1., (na+nb)))
        #%     * compute Gauss-Newton search direction
        e = (GC-g)*wf
        R = np.dot(D3.conj().T, D3)
        Vd = np.dot(D3.conj().T, e)
        if realFlag:
            R = np.real(R)
            Vd = np.real(Vd)
        
        
        gndir = linalg.solve(R, Vd)
        #%     * search along the gndir-direction *
        ll = 0.
        k = 1.
        V1 = Vcap+1.
        t1 = t
        while np.all(np.array(np.hstack((V1, ll<20.)))):
            t1 = t-np.dot(k, gndir)
            if ll == 19.:
                t1 = t
            
            
            a = np.array(np.hstack((1., t1[0:na].T)))
            b = np.array(np.hstack((np.zeros(1., nk), t1[int(na+1.)-1:na+nb].T)))
            a = apolystab(a, realFlag)
            #% Stabilizing the denominator
            t1[0:na] = a[1:na+1.].T
            GC = (np.dot(b, OM[int(indb)-1,:])/np.dot(a, OM[int(indg)-1,:])).T
            V1 = np.dot(((GC-g)*wf).conj().T, (GC-g)*wf)
            if verb:
                home
                np.disp(int2str(ll))
            
            
            k = k/2.
            ll = ll+1.
            if ll == 10.:
                gndir = np.dot(matdiv(Vd, linalg.norm(R)), length(R))
                k = 1.
            
            
            if ll == 20.:
                st = 1.
            
            
            
        if verb:
            home
            np.disp(np.array(np.hstack(('      ', getString(message('signal:invfreqs:ITERATION')), int2str(l)))))
            np.disp(np.array(np.hstack((getString(message('signal:invfreqs:CurrentFit')), num2str(V1), getString(message('signal:invfreqs:PreviousFit')), num2str(Vcap)))))
            np.disp(getString(message('signal:invfreqs:CurrentParPrevparGNdir')))
            np.disp(np.array(np.hstack((t1, t, gndir))))
            np.disp(np.array(np.hstack((getString(message('signal:invfreqs:NormOfGNvector')), num2str(linalg.norm(gndir))))))
            if st == 1.:
                np.disp(getString(message('signal:invfreqs:NoImprovement')))
                np.disp(getString(message('signal:invfreqs:IterationsThereforeTerminated')))
            
            
        
        
        t = t1
        Vcap = V1
        
    return [b, a]
def apolystab(a, realFlag):

    # Local Variables: a, realFlag, vind, v
    # Function calls: real, poly, length, apolystab, find, roots
    #%APOLYSTAB  Stabilize filter, analog
    #%   inputs: a - denominator polynomial
    #%           realFlag - 1 for real, 0 for complex
    #%   returns stabilized denoninator polynomial
    if length(a) > 0.:
        v = np.roots(a)
        vind = nonzero((np.real(v) > 0.))
        v[int(vind)-1] = -v[int(vind)-1]
        a = np.poly(v)
        if realFlag:
            a = np.real(a)
        
        
    
    
    return [a]
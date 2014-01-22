
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def arparest(x, p, method):

    # Local Variables: a, msgobj, XM, Xc, a_left, a_right, minlength_x, mx, nx, p, x, msg, e, Xc_minus, X1, method, Cz
    # Function calls: real, corrmtx, nargchk, getString, min, strcmp, issparse, nargin, length, abs, isempty, error, arparest, message, round, size
    #%ARPAREST   AR parameter estimation via a specified method.
    #%   A = ARPAREST(X,ORDER,METHOD) returns the polynomial A corresponding to 
    #%   the AR parametric signal model estimate of vector X using the specified
    #%   METHOD.  ORDER is the model order of the AR system.
    #%
    #%   Supported methods are: 'covariance' and 'modified' although all of the
    #%   methods of CORRMTX will work. In particular if 'autocorrelation' is
    #%   used, the results should be the same as those of ARYULE (but slower).
    #%
    #%   [A,E] = ARPAREST(...) returns the variance estimate E of the white noise
    #%   input to the AR model.
    #%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
    #%              Prentice-Hall, 1988, Chapter 7
    #%        S. Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATION,
    #%              Prentice-Hall, 1987, Chapter 8.
    #%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
    #%              Prentice-Hall, 1997, Chapter 3
    #%   Author(s): R. Losada and P. Pacheco
    #%   Copyright 1988-2004 The MathWorks, Inc.
    #%   $Revision: 1.5.4.3 $  $Date: 2011/05/13 18:13:56 $
    matcompat.error(nargchk(3., 3., nargin, 'struct'))
    [mx, nx] = matcompat.size(x)
    #% Initialize in case we return early
    a = np.array([])
    e = np.array([])
    #% Assign msg in case there are no errors
    msg = \'
    msgobj = np.array([])
    #% Set up necessary but not sufficient conditions for the correlation
    #% matrix to be nonsingular. From (Marple)
    _switch_val=method
    if False: # switch 
        pass
    elif _switch_val == 'covariance':
        minlength_x = 2.*p
    elif _switch_val == 'modified':
        minlength_x = 3.*p/2.
    else:
        msgobj = message('signal:arparest:UnknMethod')
        msg = getString(msgobj)
        return []
    
    #% Do some data sanity testing
    if isempty(x) or length(x)<minlength_x or matcompat.max(mx, nx) > 1.:
        if strcmp(method, 'modified'):
            msgobj = message('signal:arparest:TooSmallForModel', 'X', '3/2')
            msg = getString(msgobj)
        else:
            msgobj = message('signal:arparest:TooSmallForModel', 'X', '2')
            msg = getString(msgobj)
            
        
        return []
    
    
    if issparse(x):
        msgobj = message('signal:arparest:InputSignalCannotBeSparse')
        msg = getString(msgobj)
        return []
    
    
    if isempty(p) or p != np.round(p):
        msgobj = message('signal:arparest:ModelOrderMustBeInteger')
        msg = getString(msgobj)
        return []
    
    
    x = x.flatten(1)
    #% Generate the appropriate data matrix
    XM = corrmtx(x, p, method)
    Xc = XM[:,1:]
    X1 = XM[:,0]
    #% Coefficients estimated via the covariance method
    a_left = np.array(np.hstack((1.)))
    Xc_minus = -Xc
    a_right = linalg.solve(Xc_minus, X1)
    a = np.array(np.vstack((np.hstack((a_left)), np.hstack((a_right)))))
    #%a = [1; -Xc\X1];
    #% Estimate the input white noise variance
    Cz = np.dot(X1.conj().T, Xc)
    e = np.dot(X1.conj().T, X1)+np.dot(Cz, a[1:])
    #% Ignore the possible imaginary part due to numerical errors and force
    #% the variance estimate of the white noise to be positive
    e = np.abs(np.real(e))
    a = a.flatten(0)
    #% By convention all polynomials are row vectors
    #% [EOF] arparest.m
    return [a, e, msg, msgobj]

import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def convmtx(v, n):

    # Local Variables: cidx, c, x_left, ridx, m, n, x_right, mv, t, v, x, r, nv
    # Function calls: convmtx, length, ones, zeros, size, toeplitz
    #%CONVMTX Convolution matrix.
    #%   CONVMTX(C,N) returns the convolution matrix for vector C.
    #%   If C is a column vector and X is a column vector of length N,
    #%   then CONVMTX(C,N)*X is the same as CONV(C,X).
    #%   If R is a row vector and X is a row vector of length N,
    #%   then X*CONVMTX(R,N) is the same as CONV(R,X).
    #%
    #%   % Example:
    #%   %   Generate a simple convolution matrix.
    #%
    #%   h = [1 2 3 2 1];
    #%   convmtx(h,7)        % Convolution matrix
    #%
    #%   See also CONV.
    #%   Author(s): L. Shure, 47-88
    #%   	   T. Krauss, 3-30-93, removed dependence on toeplitz
    #%   Copyright 1988-2004 The MathWorks, Inc.
    #%   $Revision: 1.6.4.3 $  $Date: 2012/10/29 19:30:54 $
    [mv, nv] = matcompat.size(v)
    v = v.flatten(1)
    #% make v a column vector
    #t = toeplitz(np.array(np.vstack((np.hstack((v)), np.hstack((np.zeros((n-1.), 1.)))))), np.zeros(n, 1.))
    c = np.array(np.vstack((np.hstack((v)), np.hstack((np.zeros((n-1.), 1.))))))
    r = np.zeros(n, 1.)
    m = length(c)
    x_left = r[int(n)-1:2.:-1.]
    x_right = c.flatten(1)
    x = np.array(np.vstack((np.hstack((x_left)), np.hstack((x_right)))))
    #%x = [r(n:-1:2) ; c(:)];                 % build vector of user data
    cidx = np.arange(0., (m-1.)+1).conj().T
    ridx = np.arange(n, (1.)+(-1.), -1.)
    t = cidx[:,int(np.ones(n, 1.))-1]+ridx[int(np.ones(m, 1.))-1,:]
    #% Toeplitz subscripts
    t[:] = x[int(t)-1]
    #% actual data
    #% end of toeplitz code
    if mv<nv:
        t = t.T
    
    
    return [t]
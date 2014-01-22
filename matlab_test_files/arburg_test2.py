
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def arburg(x, p):

    # Local Variables: a, E, efp, k, ef, varargout, N, p, num, ebp, den, x, eb
    # Function calls: arburg, flipud, m, length, zeros, conj
    #%ARBURG   AR parameter estimation via Burg method.
    #%   A = ARBURG(X,ORDER) returns the polynomial A corresponding to the AR
    #%   parametric signal model estimate of vector X using Burg's method.
    #%   ORDER is the model order of the AR system.
    #%
    #%   [A,E] = ARBURG(...) returns the final prediction error E (the variance
    #%   estimate of the white noise input to the AR model).
    #%
    #%   [A,E,K] = ARBURG(...) returns the vector K of reflection 
    #%   coefficients (parcor coefficients).
    #%
    #%   % Example:
    #%   %   Estimate input noise variance for AR(4) model.
    #%
    #%   A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
    #%   % Generate noise standard deviations
    #%   % Seed random number generator for reproducible results
    #%   rng default;
    #%   noise_stdz=rand(50,1)+0.5;
    #%   for j=1:50
    #%       y=filter(1,A,noise_stdz(j)*randn(1024,1));
    #%       [ar_coeffs,NoiseVariance(j)]=arburg(y,4);
    #%   end
    #%   %Compare actual vs. estimated variances
    #%   plot(noise_stdz.^2,NoiseVariance,'k*');
    #%   xlabel('Input Noise Variance');
    #%   ylabel('Estimated Noise Variance');
    #%
    #%   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.
    #%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
    #%              Prentice-Hall, 1988, Chapter 7
    #%        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
    #%              Macmillan, 1988, Chapter 5
    #%   Author(s): D. Orofino and R. Losada
    #%   Copyright 1988-2009 The MathWorks, Inc.
    #%   $Revision: 1.12.4.7 $  $Date: 2012/10/29 19:30:37 $
    #% error(nargchk(2,2,nargin,'struct'))
    #% 
    #% % Check the input data type. Single precision is not supported.
    #% try
    #%     chkinputdatatype(x,p);
    #% catch ME
    #%     throwAsCaller(ME);
    #% end
    #% 
    #% validateattributes(x,{'numeric'},{'nonempty','finite','vector'},'arburg','X');
    #% validateattributes(p,{'numeric'},{'positive','integer','scalar'},'arburg','ORDER');
    #% if issparse(x),
    #%    error(message('signal:arburg:Sparse'))
    #% end
    #% if numel(x) < p+1
    #%     error(message('signal:arburg:InvalidDimension', p + 1));
    #% end
    x = x.flatten(1)
    N = length(x)
    #% Initialization
    ef = x
    eb = x
    a = 1.
    #% Initial error
    E = np.dot(x.conj().T, x)/N
    #% Preallocate 'k' for speed.
    k = np.zeros(1., p)
    #% for m=1:p
    #% Calculate the next order reflection (parcor) coefficient
    efp = ef[1:]
    ebp = eb[0:0-1.]
    num = np.dot(np.dot(-2., ebp.conj().T), efp)
    den = np.dot(efp.conj().T, efp)+np.dot(ebp.conj().T, ebp)
    k[int(m)-1] = num/den
    #% Update the forward and backward prediction errors
    ef = efp+np.dot(k[int(m)-1], ebp)
    eb = ebp+np.dot(k[int(m)-1].conj().T, efp)
    #% Update the AR coeff.
    a = np.array(np.vstack((np.hstack((a)), np.hstack((0.)))))+np.dot(k[int(m)-1], np.array(np.vstack((np.hstack((0.)), np.hstack((np.conj(np.flipud(a))))))))
    #% Update the prediction error
    E[int((m+1.))-1] = np.dot(1.-np.dot(k[int(m)-1].conj().T, k[int(m)-1]), E[int(m)-1])
    #% end
    #% a = a(:).'; % By convention all polynomials are row vectors
    #% varargout{1} = a;
    #% if nargout >= 2
    #%     varargout{2} = E(end);
    #% end
    #% if nargout >= 3
    #%     varargout{3} = k(:);
    #% end
    return [varargout]
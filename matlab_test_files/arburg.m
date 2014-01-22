function varargout = arburg( x, p)
%ARBURG   AR parameter estimation via Burg method.
%   A = ARBURG(X,ORDER) returns the polynomial A corresponding to the AR
%   parametric signal model estimate of vector X using Burg's method.
%   ORDER is the model order of the AR system.
%
%   [A,E] = ARBURG(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARBURG(...) returns the vector K of reflection 
%   coefficients (parcor coefficients).
%
%   % Example:
%   %   Estimate input noise variance for AR(4) model.
%
%   A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
%   % Generate noise standard deviations
%   % Seed random number generator for reproducible results
%   rng default;
%   noise_stdz=rand(50,1)+0.5;
%   for j=1:50
%       y=filter(1,A,noise_stdz(j)*randn(1024,1));
%       [ar_coeffs,NoiseVariance(j)]=arburg(y,4);
%   end
%   %Compare actual vs. estimated variances
%   plot(noise_stdz.^2,NoiseVariance,'k*');
%   xlabel('Input Noise Variance');
%   ylabel('Estimated Noise Variance');
%
%   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
%              Macmillan, 1988, Chapter 5

%   Author(s): D. Orofino and R. Losada
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.12.4.7 $  $Date: 2012/10/29 19:30:37 $

error(nargchk(2,2,nargin,'struct'))

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(x,p);
catch ME
    throwAsCaller(ME);
end

validateattributes(x,{'numeric'},{'nonempty','finite','vector'},'arburg','X');
validateattributes(p,{'numeric'},{'positive','integer','scalar'},'arburg','ORDER');
if issparse(x),
   error(message('signal:arburg:Sparse'))
end
if numel(x) < p+1
    error(message('signal:arburg:InvalidDimension', p + 1));
end

x  = x(:);
N  = length(x);

% Initialization
ef = x;
eb = x;
a = 1;

% Initial error
E = x'*x./N;

% Preallocate 'k' for speed.
k = zeros(1, p);

for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   efp = ef(2:end);
   ebp = eb(1:end-1);
   num = -2.*ebp'*efp;
   den = efp'*efp+ebp'*ebp;
   
   k(m) = num ./ den;
   
   % Update the forward and backward prediction errors
   ef = efp + k(m)*ebp;
   eb = ebp + k(m)'*efp;
   
   % Update the AR coeff.
   a=[a;0] + k(m)*[0;conj(flipud(a))];
   
   % Update the prediction error
   E(m+1) = (1 - k(m)'*k(m))*E(m);
end

a = a(:).'; % By convention all polynomials are row vectors
varargout{1} = a;
if nargout >= 2
    varargout{2} = E(end);
end
if nargout >= 3
    varargout{3} = k(:);
end

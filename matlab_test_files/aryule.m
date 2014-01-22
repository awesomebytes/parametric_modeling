function [a,e,k] = aryule( x, p)
%ARYULE   AR parameter estimation via Yule-Walker method.
%   A = ARYULE(X,ORDER) returns the polynomial A corresponding to the AR
%   parametric signal model estimate of vector X using the Yule-Walker
%   (autocorrelation) method.  ORDER is the model order of the AR system. 
%   This method solves the Yule-Walker equations by means of the Levinson-
%   Durbin recursion.
%
%   [A,E] = ARYULE(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARYULE(...) returns the vector K of reflection coefficients.
%   
%   % Example:
%   %   Estimate model order using decay of reflection coefficients.
%
%   rng default;
%   y=filter(1,[1 -0.75 0.5],0.2*randn(1024,1));
%
%   % Create AR(2) process
%   [ar_coeffs,NoiseVariance,reflect_coeffs]=aryule(y,10);
%
%   % Fit AR(10) model
%   stem(reflect_coeffs); axis([-0.05 10.5 -1 1]);
%   title('Reflection Coefficients by Lag'); xlabel('Lag');
%   ylabel('Reflection Coefficent');
%
%   See also PYULEAR, ARMCOV, ARBURG, ARCOV, LPC, PRONY.

%   Ref: S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
%              Macmillan, 1988, Chapter 5
%        M. Hayes, STATISTICAL DIGITAL SIGNAL PROCESSING AND MODELING, 
%              John Wiley & Sons, 1996, Chapter 8

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.6 $  $Date: 2012/10/29 19:30:38 $

error(nargchk(2,2,nargin,'struct'))

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(x,p);
catch ME
    throwAsCaller(ME);
end

[mx,nx] = size(x);
if isempty(x) || length(x) < p || min(mx,nx) > 1,
   error(message('signal:aryule:InvalidDimensions'));
elseif isempty(p) || ~(p == round(p))
   error(message('signal:aryule:MustBeInteger'))
end
if issparse(x)
   error(message('signal:aryule:Sparse'))
end

R = xcorr(x,p,'biased');
[a,e,k] = levinson(R(p+1:end),p);



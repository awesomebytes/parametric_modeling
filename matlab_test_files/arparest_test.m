function [a,e,msg,msgobj] = arparest( x, p, method)
%ARPAREST   AR parameter estimation via a specified method.
%   A = ARPAREST(X,ORDER,METHOD) returns the polynomial A corresponding to 
%   the AR parametric signal model estimate of vector X using the specified
%   METHOD.  ORDER is the model order of the AR system.
%
%   Supported methods are: 'covariance' and 'modified' although all of the
%   methods of CORRMTX will work. In particular if 'autocorrelation' is
%   used, the results should be the same as those of ARYULE (but slower).
%
%   [A,E] = ARPAREST(...) returns the variance estimate E of the white noise
%   input to the AR model.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATION,
%              Prentice-Hall, 1987, Chapter 8.
%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
%              Prentice-Hall, 1997, Chapter 3

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.5.4.3 $  $Date: 2011/05/13 18:13:56 $

error(nargchk(3,3,nargin,'struct'))
[mx,nx] = size(x);

% Initialize in case we return early
a = []; e = [];

% Assign msg in case there are no errors
msg ='';
msgobj = [];

% Set up necessary but not sufficient conditions for the correlation
% matrix to be nonsingular. From (Marple)
switch method,
case 'covariance',
   minlength_x = 2*p;
case 'modified',
   minlength_x = 3*p/2;
otherwise
   msgobj = message('signal:arparest:UnknMethod');
   msg = getString(msgobj);
   return
end

% Do some data sanity testing
if isempty(x) || length(x) < minlength_x || min(mx,nx) > 1,
    if strcmp(method, 'modified')
        msgobj = message('signal:arparest:TooSmallForModel','X','3/2');
        msg = getString(msgobj);
    else
        msgobj = message('signal:arparest:TooSmallForModel','X','2');
        msg = getString(msgobj);
    end
    return
end
if issparse(x),
   msgobj = message('signal:arparest:InputSignalCannotBeSparse');
   msg = getString(msgobj);
   return
end
if isempty(p) || p ~= round(p),
   msgobj = message('signal:arparest:ModelOrderMustBeInteger');
   msg = getString(msgobj);
   return
end

x  = x(:);

% Generate the appropriate data matrix
XM = corrmtx(x,p,method);
Xc = XM(:,2:end);
X1 = XM(:,1);

% Coefficients estimated via the covariance method
a_left = [1]
Xc_minus = -Xc
a_right = Xc_minus\X1
a = [a_left; a_right]
%a = [1; -Xc\X1];

% Estimate the input white noise variance
Cz = X1'*Xc;
e = X1'*X1 + Cz*a(2:end);

% Ignore the possible imaginary part due to numerical errors and force
% the variance estimate of the white noise to be positive
e = abs(real(e));

a = a(:).'; % By convention all polynomials are row vectors

% [EOF] arparest.m

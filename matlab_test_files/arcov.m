function [a,e] = arcov( x, p)
%ARCOV   AR parameter estimation via covariance method.
%   A = ARCOV(X,ORDER) returns the polynomial A corresponding to the AR
%   parametric signal model estimate of vector X using the Covariance method.
%   ORDER is the model order of the AR system.
%
%   [A,E] = ARCOV(...) returns the variance estimate E of the white noise
%   input to the AR model.
%
%   See also PCOV, ARMCOV, ARBURG, ARYULE, LPC, PRONY.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        P. Stoica and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
%              Prentice-Hall, 1997, Chapter 3

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.13.4.3 $  $Date: 2011/05/13 18:06:51 $

error(nargchk(2,2,nargin,'struct'));

[a,e,msg,msgobj] = arparest(x,p,'covariance');
if ~isempty(msg), error(msgobj); end

% [EOF] - arcov.m

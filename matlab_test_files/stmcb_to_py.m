function [b,a] = stmcb( x, u_in, q, p, niter, a_in )
%STMCB Compute linear model via Steiglitz-McBride iteration
%   [B,A] = stmcb(H,NB,NA) finds the coefficients of the system 
%   B(z)/A(z) with approximate impulse response H, NA poles and 
%   NB zeros.
%
%   [B,A] = stmcb(H,NB,NA,N) uses N iterations.  N defaults to 5.
%
%   [B,A] = stmcb(H,NB,NA,N,Ai) uses the vector Ai as the initial 
%   guess at the denominator coefficients.  If you don't specify Ai, 
%   STMCB uses [B,Ai] = PRONY(H,0,NA) as the initial conditions.
%
%   [B,A] = STMCB(Y,X,NB,NA,N,Ai) finds the system coefficients B and 
%   A of the system which, given X as input, has Y as output.  N and Ai
%   are again optional with default values of N = 5, [B,Ai] = PRONY(Y,0,NA).
%   Y and X must be the same length.
%
%   % Example:
%   %   Approximate the impulse response of a Butterworth filter with a 
%   %   system of lower order.
%
%   [b,a] = butter(6,0.2);              % Butterworth filter design
%   h = filter(b,a,[1 zeros(1,100)]);   % Filter data using above filter
%   freqz(b,a,128)                      % Frequency response 
%   [bb,aa] = stmcb(h,4,4);                  
%   figure; freqz(bb,aa,128)
%
%   See also PRONY, LEVINSON, LPC, ARYULE.

%   Author(s): Jim McClellan, 2-89
%   	       T. Krauss, 4-22-93, new help and options
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.6 $  $Date: 2012/10/29 19:32:10 $

error(nargchk(3,6,nargin,'struct'))

if length(u_in) == 1,
    if nargin == 3,
        niter = 5; p = q; q = u_in;
        a_in = prony(x,0,p);
    elseif nargin == 4,
        niter = p; p = q; q = u_in; 
        a_in = prony(x,0,p);
    elseif nargin == 5,
        a_in = niter; niter = p; p = q; q = u_in; 
    end
    u_in = zeros(size(x));
    u_in(1) = 1;         % make a unit impulse whose length is same as x
else
    if length(u_in)~=length(x),
        error(message('signal:stmcb:InvalidDimensions'))
    end
    if nargin < 6
       [b,a_in] = prony(x,0,p);
    end
    if nargin < 5
       niter = 5;
    end
end

a = a_in;
N = length(x);
for i=1:niter
   u = filter( 1, a, x );
   v = filter( 1, a, u_in );
   C1 = convmtx(u(:),p+1);
   C2 = convmtx(v(:),q+1);
   C1_minus = -C1;
   T_sub1 = C1_minus(1:N,:);
   T_sub2 = C2(1:N,:);
   T = [T_sub1 T_sub2];
   %T = [ C1_minus(1:N,:) C2(1:N,:) ];
   T_minus = -T;
   T_left = T(:,2:p+q+2);
   T_right = (T_minus(:,1));
   c = T_left \ T_right;
%    c = T(:,2:p+q+2)\(T_minus(:,1));   % move 1st column to RHS and do least-squares
   a = [1; c(1:p)];                % denominator coefficients
   b = c(p+1:p+q+1);               % numerator coefficients
end
a=a.';
b=b.';


function[x,sigma] = reflect(y,u,l)
%REFLECT Reflection transformation
%
%  [x,sigma] = reflect(y,u,l) reflection transformation as
%  described in Coleman and Li ??: x is reflected point, sigma is sign vector
%  corresponding to x. y is current point, u is vector of
%  upper bounds, l is vector of lower bounds.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:31:18 $

% INITIALIZATION
n = length(y);
w = zeros(n,1); 
x = zeros(n,1); 
sigma = zeros(n,1);
%
% Map y to x and compute the sign vector
% sigma according to the Coleman/Li paper ??. 
%
arg1 = (l == 0) & (u == 1); 
arg2 = (l == 0 ) & (u == inf);
arg3 = (l == -inf) & (u == 1);
arg4 = (l == -inf) & (u == inf);

w(arg1)=rem(abs(y(arg1)),2); 
x(arg1)=min(w(arg1),2-w(arg1));
sigma(arg1 & (w <= 2-w))=sign(y(arg1 & (w <= 2-w)));
sigma(arg1 & (w > 2-w))=-sign(y(arg1 & (w > 2-w)));

x(arg2) = abs(y(arg2)); 
sigma(arg2) = sign(y(arg2));

x(arg3 & (y <=1)) = y(arg3 & (y <=1));
arg5 = arg3 & (y <=1);
sigma(arg5)= 1;
arg6 = arg3 & (y > 1);
x(arg6) = 2- y(arg6);
sigma(arg6) = -1;
x(arg4) = y(arg4); 
sigma(arg4) = 1;
sigma = sigma + 1 - abs(sigma);



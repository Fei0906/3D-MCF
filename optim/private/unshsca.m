function[xx] = unshsca(x,l,u,DS)
%UNSHSCA Unshift and unscale
%
%  xx = UNSHSCA(x,l,u,DS); vector x is shifted and scaled to yield xx.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:31:29 $

arg2 = (l== -inf) & (u < inf); 
arg3 = (l> -inf) & (u == inf);
arg4 = (l > -inf) & (u < inf);
%
% UNSCALE
xx = full(DS*x);   % always full except in scalar case.
%
% UNSHIFT
xx(arg2) = xx(arg2) + u(arg2) - 1; 
xx(arg3) = xx(arg3) + l(arg3);
xx(arg4) = xx(arg4) + l(arg4);

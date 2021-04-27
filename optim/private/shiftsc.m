function[xstart,l,u,ds,DS,c] = shiftsc(xstart,l,u,typx,caller,mtxmpy,c,H,varargin)
%SHIFTSC Shift and scale
%
% [xstart,l,u,ds,DS,c] = shiftsc(xstart,l,u,typx,caller,mtxmpy,c,H) shift
% and scale vectors u,l, and xstart so that finite value of u map to
% unity and finite values of l map to zero.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:31:19 $

n = length(xstart);
ds = ones(n,1); 
DS = sparse(1:n,1:n,ds);

arg1 = (l== -inf) & (u == inf);
if nnz(arg1) == n, 
    return; 
end

arg2 = (l== -inf) & (u < inf); 
arg3 = (l> -inf) & (u == inf);
arg4 = (l > -inf) & (u < inf);

% SHIFT
vshift = zeros(n,1);
xstart(arg2) = xstart(arg2) + 1 - u(arg2);
vshift(arg2) = u(arg2) - 1;
u(arg2) = 1;

xstart(arg3) = xstart(arg3) - l(arg3);
vshift(arg3) = l(arg3);
l(arg3) = 0;

xstart(arg4) = xstart(arg4) - l(arg4);
vshift(arg4) = l(arg4);
u(arg4) = u(arg4) - l(arg4);
l(arg4) = 0;

if isequal(caller,'sqpbox')
    w = feval(mtxmpy,H,vshift,varargin{:});
else % sllsbox
    w = feval(mtxmpy,H,vshift,0,varargin{:});
end
c = c + w;

% SCALE
ds(arg1) = max(abs(typx(arg1)),ds(arg1));
ds(arg4) = abs(u(arg4));
DS = sparse(1:n,1:n,ds); 
if nargin > 5  
    c = DS*c; 
end
u = u./ds;
xstart = xstart./ds;
xstart = full(xstart); l = full(l); u = full(u); ds = full(ds);
c = full(c);

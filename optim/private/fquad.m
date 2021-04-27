function[val,g] = fquad(x,c,H,caller,mtxmpy,D,varargin)
%FQUAD Evaluate quadratic or linear least squares function.
%	val = FQUAD(x,c,H,caller,mtxmpy,D) evaluates the quadratic
%   function val = c'*x + .5*x'*D*MTX*D*x, where
%   D is a diagonal matrix and MTX is defined by the 
%   matrix multiply routine 'mtxmpy' and 'H'.
%
%   [val,g] = FQUAD(x,c,H,mtxmpy,D) also evaluates the 
%   gradient: g = D*MTX*D*x + c.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2010/10/08 17:13:44 $

if nargin < 5
   error(message('optim:fquad:NotEnoughInputs'))
end
n = length(x);
if nargin < 6 || isempty(D), 
   D = speye(n); 
end
w = full(D*x);                  % w always full unless both scalar
if isequal(caller,'sqpbox')
    ww = feval(mtxmpy,H,w,varargin{:});
else % sllsbox
    ww = feval(mtxmpy,H,w,0,varargin{:});
end
w = full(D*ww);                 % w always full unless both scalar
g = w + c; 
val = x'*((0.5*w) + c);



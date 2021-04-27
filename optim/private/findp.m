function P  = findp(A)
%FINDP  Nonsingular basis permutation matrix.
%   P = FINDP(A) locates a permutation matrix such that
%   the leading square matrix is nonsingular. If A has 
%   more columns (n) than rows(m), and A is of rank m, then
%   P is a permutation matrix such that the leading m-matrix
%   of A*P is of full rank.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/08/03 21:31:02 $

[m,n] = size(A);
if m < n, 
   A = A';
end
if ~issparse(A), 
   A=sparse(A); 
end;
pp = colamd(A);
[L,U,P] = lu(A(:,pp));



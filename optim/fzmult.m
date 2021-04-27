function [W,L,U,pcol,P]  = fzmult(A,V,transpose,L,U,pcol,P)
%

%FZMULT Multiplication with fundamental nullspace basis.
%   W = FZMULT(A,V) computes the product W of matrix Z with matrix V, that
%   is, W = Z*V, where Z is a fundamental basis for the nullspace of matrix
%   A. A must be a sparse m-by-n matrix where m < n, rank(A) = m, and
%   rank(A(1:m,1:m)) = m.  V must be p-by-q, where p = n-m. If V is sparse
%   W is sparse, else W is full.
%
%   W = FZMULT(A,V,'transpose') computes the product of the transpose of
%   the fundamental basis times V, that is, W = Z'*V. V must be p-by-q,
%   where q = n-m. FZMULT(A,V) is the same as FZMULT(A,V,[]).
%
%   [W,L,U,pcol,P]  = FZMULT(A,V) returns the sparse LU-factorization of
%   matrix A(1:m,1:m), that is, A1 = A(1:m,1:m) and P*A1(:,pcol) = L*U.
%
%   W = FZMULT(A,V,TRANSPOSE,L,U,pcol,P) uses the precomputed sparse LU
%   factorization of matrix A(1:m,1:m), that is, A1 = A(1:m,1:m) and
%   P*A1(:,pcol) = L*U.
%
%   The nullspace basis matrix Z is not formed explicitly. An implicit
%   representation is used based on the sparse LU factorization of
%   A(1:m,1:m).

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:05:52 $

% Initialization
if nargin < 2
   error(message('optim:fzmult:NotEnoughInputs')); 
end
[m,n] = size(A); 
if m >= n
   error(message('optim:fzmult:InvalidA')); 
end
if ~issparse(A)
   error(message('optim:fzmult:ANotSparse')); 
end
if nargin < 3 || isempty(transpose)
   transpose = '';     % default
end

[mm,p] = size(V);
switch transpose
case ''
   if mm ~= n-m
      error(message('optim:fzmult:SizeMismatch')); 
   end
case 'transpose'
   if mm ~= n
      error(message('optim:fzmult:SizeMismatch')); 
   end
otherwise
   error(message('optim:fzmult:InvalidTranspose'));
end

A1 = A(:,1:m);      % A1 is square 
A2 = A(:,m+1:n); 
if nargin  <  7
   pcol = colamd(A1); 
   [L,U,P] = lu(A1(:,pcol));
end

switch transpose
case ''
   % Disable the warnings about conditioning for singular and
   % nearly singular matrices
   warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
   warningstate2 = warning('off', 'MATLAB:singularMatrix');
   W = -A2*V; 
   WW = L\P*W; 
   W = U\WW; 
   WW(pcol,:) = W;
   W = [WW;V]; 
   % Restore the warning states to their original settings
   warning(warningstate1)
   warning(warningstate2)
case 'transpose'
   % Disable the warnings about conditioning for singular and
   % nearly singular matrices
   warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
   warningstate2 = warning('off', 'MATLAB:singularMatrix');
   V1 = V(1:m,:); 
   V2 = V(m+1:n,:);
   WW = -V1(pcol,:); 
   W = U'\WW; 
   WW = L'\W; 
   W = P'*WW; 
   WW = A2'*W;
   W = WW + V2; 
   % Restore the warning states to their original settings
   warning(warningstate1)
   warning(warningstate2)
otherwise
   error(message('optim:fzmult:InvalidTranspose'));
end

% Follow regular rules of sparse/full multiplication.
% A is always sparse. Thus V determines if W is sparse or full.
if ~issparse(V)
    W = full(W);
else
    W = sparse(W);
end







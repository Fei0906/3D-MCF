function [A,b]  = dep(A,tol,b)
%DEP  Remove dependent rows from the sparse system Ax = b.
%
%     [A,b]  = DEP(A,tol,b) by repeatedly applying lu-factorization
% 	   to A', where A has more cols than rows, remove numerically 
%     dependent rows of A and correspondingly adjust b. tol is a 
%     numerical tolerance used.	
%
%     WARNING: can be very expensive computationally.
%

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2010/10/08 17:13:38 $

%
% Initialization
if nargin <2, 
   tol = 1e-12; 
end
if isempty(tol); 
   tol = 1e-12; 
end
[m,n] = size(A);
if m > n 
   error(message('optim:dep:InvalidA'))
end
AA = sparse(A');
[mm,nn] = size(AA);
if nargin < 3, 
   b = zeros(m,1); 
end
if isempty(b),
   b = zeros(m,1); 
end
[mb,nb] = size(b);
if mb ~= m, 
   error(message('optim:dep:SizeMismatch')) 
end
minu = 0; 
bb = b;

% LU factorization, look for 1st dependent column
while minu <= tol
   pcol = colamd(AA);
   AA = AA(:,pcol); 
   bb = bb(pcol,1);
   [L,U,P] = lu(AA);
   u = diag(U);
   [minu,indj] = min(abs(u));
   if minu < tol
      AA(:,pcol) = AA; 
      bb(pcol,:) = bb;
      j = pcol(indj);
      A1 = AA(:,1:j-1); 
      b1 = bb(1:j-1,:);
      A2 = AA(:,j+1:nn); 
      b2 = bb(j+1:nn,:);
      AA = [A1,A2]; 
      bb = [b1;b2];
      nn = nn - 1;
   end
end
AA(:,pcol) = AA; 
bb(pcol,:) = bb;
if issparse(A)  % return A to original state with less rows
   A = AA'; 
else
   A = full(AA');
end
if nargout > 1, 
   b = bb; 
end

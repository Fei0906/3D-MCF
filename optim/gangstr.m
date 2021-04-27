function A = gangstr(M,tol)
%

%GANGSTR Zero out 'small' entries subject to structural rank.
%   A = GANGSTR(M,TOL) creates matrix A of full structural rank such that A
%   is M except that elements of M that are relatively 'small', based on
%   TOL, are zeros in A. TOL may be decreased, if needed, until SPRANK(A) =
%   SPRANK(M). M must have at least as many columns as rows. Default TOL is
%   1e-2.
%
%   Elements of M relatively less than TOL are identified by first
%   normalizing all the rows of M to have norm 1. Then nonzeros in M are
%   examined in a columnwise fashion, removing the ones with values of
%   magnitude less than TOL*(maximum absolute value in that column).
%
%   See also SPRANK, SPY.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:05:53 $


if isempty(M)
   A=M;
   return
end
[m,n] = size(M);
if m > n
   error(message('optim:gangstr:InvalidMatrix'))
end
if nargin < 2, 
   tol = 1e-2; 
end

% Normalize M
Msqr=M.*M;
% Sum rows of Msqr and return column vector.
% Msqr cannot be a column vector.
X = sum(Msqr,2); % spdiags needs a column vector as input

X(X==0) = 1;
normM = spdiags(1./sqrt(X),0,m,m)*M;
A_1 = normM(:,1:m); 
A_2 = normM(:,m+1:n);

% Remove nonzeros from A_2 first
[I_2,J_2,V_2] = find(A_2);
if m > 1
    maxvec = full(tol*max(abs(A_2))); 
    maxvec = maxvec(:); 
else 
    maxvec = full(tol*abs(A_2)); 
end
tobekept = find(abs(V_2) > maxvec(J_2));
A_2 = sparse(I_2(tobekept),J_2(tobekept),V_2(tobekept),m,n-m);

% Remove nonzeros from A_1 making sure it has full structural rank
dim = sprank(A_1);
sprA = 0;
[I_1,J_1,V_1] = find(A_1);
absA_1 = abs(A_1);
while sprA < dim
   % Take the max of the columns of abs(M) unless M is one row only   
   if m > 1 
      maxvec = full(tol*max(absA_1)); 
      maxvec = maxvec(:);    % make maxvec consistent with I,J,V
   else 
      maxvec = full(tol*absA_1); % maxvec, I,J,V are all rows in this case 
   end
   tobekept = find(abs(V_1) > maxvec(J_1));
   A_1 = sparse(I_1(tobekept),J_1(tobekept),V_1(tobekept),m,m);
   sprA = sprank(A_1);
   tol = tol-0.025;
end

% Denormalize
A = spdiags(sqrt(X),0,m,m)*[A_1 A_2];
if ~issparse(M)
    A = full(A);
end


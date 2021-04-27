function [Xsub,idx]=licols(X,tol)
%=======================================================================
%   Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%   Input:
%       X: The given input matrix
%       tol: A rank estimation tolerance. Default=1e-10
%
%   Output:
%       Xsub: The extracted columns of X
%       idx:  The indices (into X) of the extracted columns
%=======================================================================

%   X has no non-zeros and hence no independent columns
if ~nnz(X) 
	Xsub=[]; idx=[];
	return
end

if nargin<2, 
    tol=1e-10; 
end

[~, R, E] = qr(X,0); 
if ~isvector(R)
	diagr = abs(diag(R));
else
	diagr = R(1);   
end

%	Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last');
idx=sort(E(1:r));
Xsub=X(:,idx);

end
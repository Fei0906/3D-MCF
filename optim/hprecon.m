function[R,pvec] = hprecon(H,upperbandw,DM,DG,varargin)
%

%HPRECON Sparse Cholesky factor of H-preconditioner
%
%   [R,PVEC] = HPRECON(H,UPPERBANDW,DM,DG) computes the 
%   sparse Cholesky factor (transpose of a (usually) banded 
%   preconditioner of square matrix M 
%                M = DM*H*DM + DG
%   where DM and DG are non-negative sparse diagonal matrices. 
%   R'*R approximates M(pvec,pvec), i.e.
%          R'*R = M(pvec,pvec)
%
%   H may not be the true Hessian.  If H is the same size as the
%   true Hessian, H will be used in computing the preconditioner R.
%   Otherwise, compute a diagonal preconditioner for
%               M = DM*DM + DG
%
%   If 0 < UPPERBANDW <  n then the upper bandwidth of 
%   R is UPPERBANDW. If UPPERBANDW >= n then the structure of R 
%   corresponds to a sparse Cholesky factorization of H
%   using the symamd ordering (the ordering is returned in PVEC).
%

%   Default preconditioner for SFMINBX and SQPMIN.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:05:57 $

if nargin < 1, 
   error(message('optim:hprecon:NotEnoughInputs')); 
end
if nargin <2, 
    upperbandw = 0; 
    if nargin < 3
        DM = [];
        if nargin < 4
            DG = [];
        end, end, end

[rows,cols] = size(H); 
n = length(DM);
% In case "H" isn't really H, but something else to use with HessMult function.
if ~isnumeric(H) || ~isequal(n,rows) || ~isequal(n,cols)
    % H is not the right size; ignore requested bandwidth and compute
    % diagonal preconditioner based only on DM and DG.
    pvec = (1:n);
    d1 = full(diag(DM));  % full vector
    d2 = full(diag(DG)); 
    dd = sqrt(d1.*d1 + abs(d2));
    R = sparse(1:n,1:n,dd);
    return
end

H = DM*H*DM + DG;
pvec = (1:n);
epsi = .0001*ones(n,1);
info = 1;

if upperbandw >= n-1 % Try complete approximation to H
   pvec = symamd(H);
   ddiag = diag(H);
   mind = min(ddiag);
   lambda = 0;
   if mind < 0, 
      lambda = -mind + .001; 
   end
   while info > 0
      H = H + lambda*speye(n);
      [R,info] = chol(H(pvec,pvec));
      lambda = lambda + 10;
   end
elseif (upperbandw > 0) && ( upperbandw < n-1) % Banded approximation to H
   % Banded approximation
   lambda = 0;
   ddiag = diag(H);
   mind = min(ddiag);
   if mind < 0, 
      lambda = -mind + .001; 
   end
   H = tril(triu(H,-upperbandw),upperbandw);
   while info > 0
      H = H + lambda*speye(n);
      [R,info] = chol(H);
      lambda = 4*lambda;
      if lambda <= .001, 
         lambda = 1; 
      end
   end
elseif upperbandw == 0 % diagonal approximation for H
   dnrms = sqrt(sum(H.*H))';
   d = max(sqrt(dnrms),epsi);
   R = sparse(1:n,1:n,full(d));
   pvec = (1:n);
else
    error(message('optim:hprecon:InvalidUpperbandw'))
end


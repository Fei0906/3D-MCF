function[RPCMTX,ppvec] = pceye(A,pcoptions,DM,DG,varargin)
%

%PCEYE Precondition based on DM and DG.
% Produce diagonal preconditioner (factor) for
%
%       M = DM*(A'*A)*DM + DG
%
% where DM and DG are non-negative sparse diagonal matrices, and A' is
% unknown (empty input argument).
%

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/09 01:06:06 $

  if nargin < 2
    error(message('optim:pceye:NotEnoughInputs'))
  end
  n = length(DM);
  ppvec = (1:n)';
  d1 = full(diag(DM)); 
  d2 = full(diag(DG)); 
  dd = sqrt(d1.*d1 + abs(d2));
  RPCMTX = sparse(1:n,1:n,dd);

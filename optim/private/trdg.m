function[s,snod,qpval,pdef,pcgit,Z]=trdg(x,g,H,...
   delt,dv,mtxmpy,pcmtx,pcf,ptol,kmax,A,b,Z,dnewt,options,defaultopt,...
   PT,LZ,UZ,pcolZ,PZ,varargin)
%TRDG Trust region dogleg step for linear equality problem
%
% [s,snod,qpval,pdef,pcgit,RP,pp,R,p,Z]=trdg(x,g,H,...
%	       delt,dv,mtxmpy,pcmtx,pcf,ptol,th,A,b,Z,dnewt) determines
% a trust region dogleg step for linear equality problem. The
% driver function is SFMINLE
%
% Determine the trial step `s', an approx. trust region solution.
% `s' is chosen as the best of 3 steps: the scaled gradient
%  (truncated to  maintain strict feasibility),
%  a 2-D trust region solution (truncated to remain strictly feas.),
%  and the reflection of the 2-D trust region solution,
%  (truncated to remain strictly feasible).

%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/04/15 23:21:32 $

% Initialization
[m,n] = size(A);  tol = sqrt(eps);
pcgit = 0;  nrmdg = norm(g); D = speye(n);

DG = sparse(n,n);
if nargin < 12, ptol = .1; end
v1 = dnewt; pdef = 1;

% Determine 2-D subspace Z?
if isempty(Z); 
   if isempty(v1);
      %kmax = floor(n/2);
      verb = 0; % Do not print messages we need from ppcgr when called from sqpmin
      computeLambda = 0;  % Do not need lagrange multipliers here
      [v1,pdef,pcgit]=ppcgr(g,H,mtxmpy,A,b,options,defaultopt,verb,computeLambda,ptol,PT,LZ,UZ,pcolZ,PZ,varargin{:});
   end
   v1 = v1/norm(v1); Z(:,1) = v1;
   if (pdef < 1) % Not positive definite
      v2 = D*sign(g); v2 = v2/norm(v2);
      v2 = project(A,v2,PT,LZ,UZ,pcolZ,PZ);
      v2 = v2 - v1*(v1'*v2); nrmv2 = norm(v2);
      if nrmv2 > tol, v2 = v2/nrmv2; Z(:,2) = v2; end
   else % Positive definite
      v2 = g/norm(g);
      v2 = project(A,v2,PT,LZ,UZ,pcolZ,PZ);
      v2 = v2 - v1*(v1'*v2); nrmv2 = norm(v2);
      if nrmv2 > tol, v2 = v2/nrmv2; Z(:,2) = v2; end
   end
end

% Reduce to the 2-D subspace Z
W = D*Z; 
WW = feval(mtxmpy,H,W(1:n,:),varargin{:});
[mWW,nWW] = size(WW);
W = D*WW;
MM = full(Z'*W + Z'*DG*Z); rhs=full(Z'*g);
[st,qpval,po,fcnt,lam] = trust(rhs,MM,delt);
ss = Z*st;
s = abs(diag(D)).*ss;
s = full(s);  snod = s;
qpval1 = rhs'*st + (.5*st)'*MM*st;



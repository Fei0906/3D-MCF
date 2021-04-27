function[s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,...
  pcmtx,pcoptions,tol,H,llsprob,kmax,varargin)
%DRQPBOX Descent direction for quadratic problem.
%	[s,posdef,pcgit] = DRQPBOX(D,DS,grad,delta,g,dv,mtxmpy,...
%                                   pcmtx,pcoptions,tol,H,llsprob)
%   determines s, a descent direction (for use with SLLSBOX,SQPBOX) 
%   for quadratic minimization subject to box constraints.
%   If negative curvature is discovered in the CG process
%   then posdef = 0; otherwise posdef = 1. pcgit is the
%   number of CG iterations. LLSPROB is a flag that tells if
%   the caller is SLLSBOX or SQPBOX so the preconditioner can be
%   called correctly.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/08/03 21:30:58 $

n = length(grad); 
tsize = 0; 
pcgit = 0; 
tol2 = sqrt(eps);

DM = D*DS;
tau = 1e-4; 
vtau = tau*ones(n,1); 
diagDM = full(diag(DM));
ddiag = abs(g).*dv; 
arg = (abs(g) + diagDM < tau);
ddiag(arg == 1) = vtau(arg == 1);
DG = sparse(1:n,1:n,full(ddiag));

% A PRECONDITIONED CONJUGATE GRADIENT ROUTINE IS USED.
if llsprob
   [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
   [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,mtxmpy,H,R,permR,'jacobprecon',pcoptions,varargin{:});
else
   % preconditioner based on H, no matter what it is
   [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
   [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,mtxmpy,H,R,permR,'hessprecon',pcoptions,varargin{:});
end

% FORM A 2-D SUBSPACE
if norm(v1) == 0
   % Special treatment if CG step is zero
   s = zeros(n,1);
else
   v1 = v1/norm(v1);
   Z(:,1) = v1;
   if (posdef < 1)
      v2 = D*sign(grad);
      v2 = v2/norm(v2);
      v2 = v2 - v1*(v1'*v2);
      nrmv2 = norm(v2);
      if nrmv2 > tol2
         v2 = v2/nrmv2; Z(:,2) = v2;
      end
   else
      v1 = v1/norm(v1); 
      v2 = grad/norm(grad);
      v2 = v2 - v1*(v1'*v2); 
      nrmv2 = norm(v2);
      if nrmv2 > tol2
         v2 = v2/nrmv2; Z(:,2) = v2;
      end
   end

   W = DM*Z; 
   if llsprob
      WW = feval(mtxmpy,H,W,0,varargin{:}); 
   else
      WW = feval(mtxmpy,H,W,varargin{:}); 
   end
   W = DM*WW;
   MM = Z'*W + Z'*DG*Z;
   rhs=full(Z'*grad);

   % SOLVE TRUST REGION OVER 2-D.
   [ss,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
   ss = Z*ss;
   norms = norm(ss);
   s = abs(diag(D)).*ss;
end





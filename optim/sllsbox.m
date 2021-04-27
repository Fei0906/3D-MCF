function [x,residnorm,residual,fnrm,it,npcg,exitflag,LAMBDA,msg] = ...
    sllsbox(A,b,lb,ub,xstart,params,options,defaultopt,mtxmpy,varargin)
%

%SLLSBOX Linear least-squares with bounds
%
% x=slls(A,b,lb,ub) returns the solution to the
% box-constrained linear least-squares problem,
%
%        min { ||Ax - b||^2 : lb <= x <= ub}.
%
% where A is a matrix with more rows than cols. (may be virtual)
%

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:06:14 $

if nargin < 4
  error(message('optim:sllsbox:NotEnoughInputs'))
end
if nargin < 5
    xstart = []; 
end
if nargin < 6 
    params.verb = 1;
    params.emptyInitialPoint = true;
end 
if nargin < 7
    options = []; 
end

pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tol = optimget(options,'TolFun',defaultopt,'fast') ;
itb = optimget(options,'MaxIter',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast');

pcmtx = optimget(options,'Preconditioner',@aprecon) ; % no default, use slow optimget

c = feval(mtxmpy,A,-b,-1,varargin{:}); 
n = length(c); 
it = 1; 
cvec = c; 
% In case the defaults were gathered from calling: optimset('lsqlin'):
numberOfVariables = n;
kmax = optimget(options,'MaxPCGIter',defaultopt,'fast') ;
typx = optimget(options,'TypicalX',defaultopt,'fast') ;
if ischar(kmax)
   if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
      kmax = max(1,floor(numberOfVariables/2));
   else
      error(message('optim:sllsbox:InvalidMaxPCGIter'))
   end
end
if ischar(typx)
   if isequal(lower(typx),'ones(numberofvariables,1)')
      typx = ones(numberOfVariables,1);
   else
      error(message('optim:sllsbox:InvalidTypicalX'))
   end
end
checkoptionsize('TypicalX', size(typx), numberOfVariables);

if n == 0
   error(message('optim:sllsbox:InvalidN'))
end
if isempty(lb), lb = -inf*ones(n,1); end
if isempty(ub), ub = inf*ones(n,1); end
arg = (ub >= 1e10); arg2 = (lb <= -1e10);
ub(arg) = inf;
lb(arg2) = -inf;
if any(ub == lb) 
   error(message('optim:sllsbox:EqualLowerUpperBnd'))
elseif min(ub-lb) <= 0
   error(message('optim:sllsbox:InconsistentBnds'))
end
lvec = lb; uvec = ub;

% Initial point
xinitOutOfBounds_idx = xstart < lb | xstart > ub;
if any(xinitOutOfBounds_idx)
    if params.emptyInitialPoint
        % When no x0 provided, x0 is set to the zero vector in lsqlin.
        % If zero start point violates at least one bound, set x0
        % to box-centered point. (This is done this way to maintain
        % backwards compatibility.)
        xstart = startx(ub,lb);
    else
        % If user-provided x0 has components not within bounds,
        % set those components to a box-centered point
        xstart = startx(ub,lb,xstart,xinitOutOfBounds_idx);
    end
end

if isempty(typx), typx = ones(n,1); end
if isempty(params.verb), params.verb = 1; end

%   SHIFT AND SCALE
[xstart,lb,ub,ds,DS,c] = shiftsc(xstart,lb,ub,typx,'sllsbox',mtxmpy,cvec,A,varargin{:});

%   MORE INITIALIZATIONS
tol1 = tol; tol2 = sqrt(tol1);
dellow = 1.; delup = 10^3; npcg = 0; digits = inf; 
done = false;
del = 10*eps; posdef = 0;x = xstart; y = x;sigma = ones(n,1);
oval = inf; [val,g] = fquad(x,c,A,'sllsbox',mtxmpy,DS,varargin{:});
prev_diff = 0;

%   MAIN LOOP: GENERATE FEAS. SEQ.  x(it) S.T. q(x(it)) IS DECREASING.
while ~done
   %     Update and display
   [v,dv] = definev(g,x,lb,ub);
   fnrm = norm(v.*g,inf); 
   %
   %     TEST FOR CONVERGENCE
   diff = abs(oval-val);
   if it > 1, digits = (prev_diff)/max(diff,eps); end
   prev_diff = diff; oval = val; 
   if diff < tol1*(1+abs(oval)),
      exitflag = 3; done = true;
      msg = sprintf(['Optimization terminated: relative function value changing by less\n' ...
                     ' than OPTIONS.TolFun.']);
      if params.verb > 0
         disp(msg)
      end
   elseif (  (diff < tol2*(1+abs(oval))) && (digits < 3.5)),
      exitflag = 3; done = true;
      msg = sprintf(['Optimization terminated: relative function value changing by less\n' ...
                     ' than sqrt(OPTIONS.TolFun), and rate of progress is slow.']);
      if params.verb > 0
         disp(msg)
      end
      
   elseif ((fnrm < tol1) && (posdef ==1)) 
      exitflag = 1; done = true;
      msg = sprintf(['Optimization terminated: first order optimality with optimality\n' ...
                     ' gradient norm less than OPTIONS.TolFun.']);
      if params.verb > 0
         disp(msg)
      end
      
   end
   %
   if ~done
      %       DETERMINE THE SEARCH DIRECTION
      dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd).*sigma));
      grad = D*g; normg = norm(grad); 
      delta = max(dellow,norm(v)); delta = min(delta,delup);
      [s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,...
         pcmtx,pcflags,pcgtol,A,1,kmax,varargin{:});
      npcg = npcg + pcgit;
      
      %       DO A REFLECTIVE (BISECTION) LINE SEARCH. UPDATE x,y,sigma.
      strg= s'*(sigma.*g); ox = x;  osig = sigma; ostrg = strg;
      if strg >= 0, 
         exitflag = -4; done = true;
         msg = sprintf(['Optimization terminated: ill-conditioning prevents' ...
                        ' further optimization.']);
         if params.verb > 0
            disp(msg)
         end
      else
         [x,sigma,alpha] = biqpbox(s,c,ostrg,ox,y,osig,lb,ub,oval,posdef,...
            normg,DS,mtxmpy,A,1,varargin{:});
         y = y + alpha*s; 
         
         %          PERTURB x AND y ?
         [pert,x,y] = perturbTrustRegionReflective(x,lb,ub,del,y,sigma);
         
         %          EVALUATE NEW FUNCTION VALUE, GRADIENT. 
         it=it+1; [val,g] = fquad(x,c,A,'sllsbox',mtxmpy,DS,varargin{:}); 
      end
   end
   if it >= itb, 
      it=it -1; 
      exitflag = 0; done = true;
      msg = sprintf('Maximum number of iterations exceeded; increase OPTIONS.MaxIter.');
      if params.verb > 0
         disp(msg)
      end     
   end
end

%   RESCALE, UNSHIFT, AND EXIT.
x = unshsca(x,lvec,uvec,DS);

residual = b-feval(mtxmpy,A,x,1,varargin{:}); % b-A*x
residnorm = sum(residual.*residual);

g = feval(mtxmpy,A,-residual,-1,varargin{:}); % A'*(A*x-b) so -residual
g = full(g);
LAMBDA.lower = zeros(length(lvec),1);
LAMBDA.upper = zeros(length(uvec),1);
active_tol = sqrt(eps);
argl = logical(abs(x-lvec) < active_tol);
argu = logical(abs(x-uvec) < active_tol);

LAMBDA.lower(argl) = (g(argl));
LAMBDA.upper(argu) = -(g(argu));
LAMBDA.ineqlin = []; 
LAMBDA.eqlin = [];

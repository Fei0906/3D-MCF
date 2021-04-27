function [x,val,csnrm,it,npcg,exitflag,LAMBDA,msg] = sqpbox(c,H,mtxmpy,l,u,xstart,options,defaultopt,...
    numberOfVariables,verb,computeLambda,varargin)
%SQPBOX Minimize box-constrained quadratic function
%
%   Locate a (local) solution to the box-constrained QP:
%
%        min { q(x) = .5x'Hx + c'x :  l <= x <= u}.
%
%   where H is sparse symmetric, c is a col vector,
%   l,u are vectors of lower and upper bounds respectively.

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2011/10/15 01:57:48 $

% Retrieving options

% Add Preconditioner to defaultopt. This is not a documented option for quadprog, 
% and is not added to defaultopt at the user-facing function level.
defaultopt.Preconditioner = @hprecon;

pcmtx = optimget(options,'Preconditioner',defaultopt,'fast');
kmax = optimget(options,'MaxPCGIter',defaultopt,'fast');
typx = optimget(options,'TypicalX',defaultopt,'fast');
if ischar(kmax)
    if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
        kmax = max(1,floor(numberOfVariables/2));
    elseif isequal(lower(kmax),'numberofvariables')
        kmax = numberOfVariables;
    else
        error(message('optim:sqpbox:InvalidMaxPCGIter'))
    end
end
if ischar(typx)
    if isequal(lower(typx),'ones(numberofvariables,1)')
        typx = ones(numberOfVariables,1);
    else
        error(message('optim:sqpbox:InvalidTypicalX'))
    end
end
checkoptionsize('TypicalX', size(typx), numberOfVariables);
pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tolx = optimget(options,'TolX',defaultopt,'fast') ;
tolfun = optimget(options,'TolFun',defaultopt,'fast');
pcgtol = optimget(options,'TolPCG',defaultopt,'fast') ;  % pcgtol = .1;
itb = optimget(options,'MaxIter',defaultopt,'fast') ;

%   INITIALIZATIONS
if nargin <= 1
    error(message('optim:sqpbox:NotEnoughInputs'))
end
n = length(c);
it = 0;
cvec = c;

if n == 0
    error(message('optim:sqpbox:InvalidN'))
end
if nargin <= 2
    l = -inf*ones(n,1);
end
if nargin <= 3,
    u = inf*ones(n,1);
end
if isempty(l),
    l = -inf*ones(n,1);
end
if isempty(u),
    u = inf*ones(n,1);
end
arg = (u >= 1e10);
arg2 = (l <= -1e10);
u(arg) = inf;
l(arg2) = -inf;
if min(u-l) <= 0
    error(message('optim:sqpbox:InconsistentBnds'))
end
lvec = l; uvec = u;

pcgit = 0;
tolfun2 = sqrt(tolfun);
[xstart,l,u,ds,DS,c] = shiftsc(xstart,l,u,typx,'sqpbox',mtxmpy,cvec,H,varargin{:});
dellow = 1.;
delup = 10^3;
npcg = 0;
digits = inf;
done = false;
v = zeros(n,1);
dv = ones(n,1);
del = 10*eps;
posdef = 1;
x = xstart;
y = x;
sigma = ones(n,1);
g = zeros(n,1);
oval = inf;
prev_diff = 0;
[val,g] = fquad(x,c,H,'sqpbox',mtxmpy,DS,varargin{:});
[v,dv] = definev(g,x,l,u);
csnrm = norm(v.*g,inf);
if csnrm == 0
    % If initial point is a 1st order point then randomly perturb the
    % initial point a little while keeping it feasible and reinitialize.
    dir = zeros(n,1);
    pos = u-x > x-l;
    neg = u-x <= x-l;
    dir(pos) = 1; dir(neg) = -1;
    % Get random noise, but "put it back" so we don't affect anyone
    dflt = RandStream.getGlobalStream;
    randstate = dflt.State;
    x = x + dir.*rand(n,1).*max(u-x,x-l).*1e-1;
    dflt.State = randstate;
    y = x;
    [val,g] = fquad(x,c,H,'sqpbox',mtxmpy,DS,varargin{:});
end

%
%   MAIN LOOP: GENERATE FEAS. SEQ.  x(it) S.T. q(x(it)) IS DECREASING.
while ~done
    it = it + 1;
    %     Update and display
    [v,dv] = definev(g,x,l,u);
    csnrm = norm(v.*g,inf);
    r = abs(min(u-x,x-l));
    
    delta = max(dellow,norm(v));
    delta = min(delta,delup);

    %
    %     TEST FOR CONVERGENCE
    diff = abs(oval-val);
    if it > 1,
        digits = (prev_diff)/max(diff,eps);
    end
    prev_diff = diff;
    oval = val;
    if diff < tolfun*(1+abs(oval)),
        exitflag = 3; done = true;
        msg = sprintf(['Optimization terminated: relative function value changing by\n' ...
            ' less than OPTIONS.TolFun.']);
        if verb > 0
            disp(msg)
        end
    elseif ((diff < tolfun2*(1+abs(oval))) & (digits < 3.5)) & posdef,
        exitflag = 3; done = true;
        msg = sprintf(['Optimization terminated: relative function value changing by less\n' ...
            ' than sqrt(OPTIONS.TolFun), no negative curvature detected in current\n' ...
            ' trust region model and the rate of progress (change in f(x)) is slow.']);
        if verb > 0
            disp(msg)
        end
    elseif ((csnrm < tolfun) & posdef & it > 1),
        exitflag = 1; done = true;
        msg = sprintf(['Optimization terminated: no negative curvature detected in current\n' ...
            ' trust region model and first order optimality measure < OPTIONS.TolFun.']);
        if verb > 0
            disp(msg)
        end
    end
    %
    if ~done            
        %       DETERMINE THE SEARCH DIRECTION
        dd = abs(v);
        D = sparse(1:n,1:n,full(sqrt(dd).*sigma));
        grad = D*g;
        normg = norm(grad);     

        [s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,...
            pcmtx,pcflags,pcgtol,H,0,kmax,varargin{:});

        npcg = npcg + pcgit;
        %
        %       DO A REFLECTIVE (BISECTION) LINE SEARCH. UPDATE x,y,sigma.
        strg= s'*(sigma.*g);
        ox = x;
        osig = sigma;
        ostrg = strg;
        if strg >= 0,
            exitflag = -2; done = true;
            msg = sprintf(['Optimization terminated: loss of feasibility with respect to the\n' ...
                ' constraints detected.']);
            if verb > 0
                disp(msg);
            end
        else
            [x,sigma,alpha] = biqpbox(s,c,ostrg,ox,y,osig,l,u,oval,posdef,...
                normg,DS,mtxmpy,H,0,varargin{:});
            if alpha == 0,
                exitflag = -4; done = true;
                msg = sprintf(['Optimization terminated: current direction not descent direction;\n' ...
                    ' the problem may be ill-conditioned.']);
                if verb > 0
                    disp(msg)
                end

            end
            y = y + alpha*s;
            %
            %          PERTURB x AND y ?
            [~,x,y] = perturbTrustRegionReflective(x,l,u,del,y,sigma);
            %
            %          EVALUATE NEW FUNCTION VALUE, GRADIENT.
            [val,g] = fquad(x,c,H,'sqpbox',mtxmpy,DS,varargin{:});
        end
        if it >= itb,
            exitflag = 0; done = true;
            msg = sprintf('Maximum number of iterations exceeded; increase options.MaxIter');
            if verb > 0
                disp(msg)
            end
        end
    end
end

%
%   RESCALE, UNSHIFT, AND EXIT.
x = unshsca(x,lvec,uvec,DS);
% unscaled so leave out DS
[val,g] = fquad(x,cvec,H,'sqpbox',mtxmpy,[],varargin{:});

if computeLambda
    LAMBDA.lower = zeros(length(lvec),1);
    LAMBDA.upper = zeros(length(uvec),1);
    active_tol = sqrt(eps);
    argl = logical(abs(x-lvec) < active_tol);
    argu = logical(abs(x-uvec) < active_tol);

    g = full(g);
    LAMBDA.lower(argl) = abs(g(argl));
    LAMBDA.upper(argu) = -abs(g(argu));
else
    LAMBDA=[];
end


function kktError = computeKKTErrorForQPLP(H,f,Ain,bin,Aeq,beq,lb,ub,lambdas,x)
%COMPUTEKKTERRORFORQPLP computes first order optimality of a QP or an LP.
%   This function is private to linprog and quadprog.

%   Copyright 2010 The MathWorks, Inc.

% Return if the decision variable is empty.
if isempty(x)
    kktError = [];
    return;
end

nvars = length(x);
% Make sure f and x are columns.
f = f(:);
x = x(:);

if isempty(Ain), Ain = zeros(0,nvars); end
if isempty(Aeq), Aeq = zeros(0,nvars); end
mIneq = size(Ain,1); 
mEq = size(Aeq,1);

% Make sure that lambdas structure is not empty or else populate it.
if isempty(lambdas.lower)
    lambdas.lower = zeros(nvars,1);
end
if isempty(lambdas.upper)
    lambdas.upper = zeros(nvars,1);
end
if isempty(lambdas.ineqlin)
    lambdas.ineqlin = zeros(mIneq,1);
end
if isempty(lambdas.eqlin)
    lambdas.eqlin = zeros(mEq,1);
end

% Calculate the gradient of the objective.
if isempty(H)
   grad = f; 
else
   grad = H*x + f;
end
% Calculate the norm of the grad(Lagrangian) where grad(Lagrangian) is:
% (grad + Aeq'*Lambda_eq + Ain'Lambda_ineq + Alb'*Lambda_lb + Aub'*Lambda_ub)
normgradLag = norm(grad + ...
                   Aeq'*lambdas.eqlin + ...
                   Ain'*lambdas.ineqlin + ...
                  -lambdas.lower + ... % -eye(nvars)*lambdas.lower
                   lambdas.upper, Inf); % eye(nvars)*lambdas.upper

% Now compute the complementarity (c_i*Lambda_i) for all the inequalities

% Find the indices of finite bounds
finiteLB = ~isinf(lb); lenLB = nnz(finiteLB);
finiteUB = ~isinf(ub); lenUB = nnz(finiteUB);
% Initialization
Comp = zeros(mIneq+lenLB+lenUB,1);

% Inequality constraints
if mIneq > 0
   Comp(1:mIneq) = lambdas.ineqlin.*(Ain*x-bin);
end

% Lower bounds
if lenLB > 0
    Comp(mIneq+1:mIneq+lenLB) = lambdas.lower(finiteLB).*(lb(finiteLB) - x(finiteLB));
end
% Upper bounds
if lenUB > 0
    Comp(mIneq+lenLB+1:end) = lambdas.upper(finiteUB).*(x(finiteUB) - ub(finiteUB));
end
% Take the norm of the complementarity       
normComp = norm(Comp,Inf);

% KKT error
kktError = max(normgradLag, normComp);


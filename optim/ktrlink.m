function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = ktrlink(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,ktrOptsFile)
%KTRLINK Interface to third-party solver KNITRO (R).
%   KTRLINK uses the third-party solver KNITRO (R) which is
%   formulated for problems of the form:
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%     X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                              LB <= X <= UB        (bounds)
%                                  
%   KNITRO (R) is available from Ziena Optimization, Inc. For more
%   information, please consult:
%       http://www.ziena.com/
%
%   Prior to executing KTRLINK, you must first install KNITRO (R) on
%   your system. For more information on setup, see the KTRLINK
%   documentation.
%
%   X = KTRLINK(FUN,X0,A,B) starts at X0 and finds a minimum X to the
%   function FUN, subject to the linear inequalities A*X <= B. FUN accepts
%   input X and returns a scalar function value F evaluated at X. X0 may be
%   a scalar, vector, or matrix. 
%
%   X = KTRLINK(FUN,X0,A,B,Aeq,Beq) minimizes FUN subject to the 
%   linear equalities Aeq*X = Beq as well as A*X <= B. (Set A = [] and 
%   B = [] if no inequalities exist.)
%
%   X = KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in 
%   the range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. KTRLINK
%   minimizes FUN such that C(X) <= 0 and Ceq(X) = 0. (Set LB = [] and/or 
%   UB = [] if no bounds exist.)
%
%   X = KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes
%   with the default optimization parameters replaced by values in the
%   structure OPTIONS, an argument created with the OPTIMSET function. See
%   OPTIMSET for details. For a list of options accepted by FMINCON refer
%   to the documentation. Use OPTIONS = [] as a  place holder if no options
%   are set.
%  
%   X = KTRLINK(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS,
%   knitroOptionsFile) minimizes with the default optimization 
%   parameters replaced by the values in the KNITRO (R) options text file 
%   whose name is given in the string knitroOptionsFile.
%
%   [X,FVAL] = KTRLINK(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = KTRLINK(FUN,X0,...) returns an EXITFLAG that
%   describes the exit condition of KTRLINK. For possible values of
%   EXITFLAG and the corresponding exit conditions, refer to KNITRO (R) 
%   documentation.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = KTRLINK(FUN,X0,...) returns a 
%   structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm, the 
%   first-order optimality in OUTPUT.firstorderopt, and the exit message in
%   OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = KTRLINK(FUN,X0,...) returns the
%   Lagrange multipliers at the solution X: LAMBDA.lower for LB,
%   LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear inequalities,
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear
%   equalities.
%
%   KNITRO (R) is a third-party library solver available from Ziena
%   Optimization, Inc. For more information, please consult:
%       http://www.ziena.com/
% 
%   See also FMINCON, @, FUNCTION_HANDLE.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.10.19 $  $Date: 2011/08/13 17:31:19 $

defaultopt = struct( ...
    'Algorithm','interior-point', ...
    'AlwaysHonorConstraints','bounds', ...
    'Display','final', ...
    'FinDiffType','forward', ...
    'FunValCheck','off', ... 
    'GradConstr','off',...
    'GradObj','off', ...
    'HessFcn',[],...   
    'Hessian','bfgs', ...   % BFGS Hessian by default
    'HessMult',[], ...      % HessMult [] by default
    'HessPattern','sparse(ones(numberOfVariables))',...
    'InitBarrierParam',0.1, ...
    'InitTrustRegionRadius','sqrt(numberOfVariables)',...
    'JacobPattern','sparse(ones(Jrows,Jcols))', ... % Full Jacobian by default
    'MaxIter',10000, ...
    'MaxProjCGIter','2*(numberOfVariables-numberOfEqualities)', ...
    'ObjectiveLimit',1e20, ...
    'ScaleProblem','obj-and-constr', ...
    'SubproblemAlgorithm','ldl-factorization',...
    'TolCon',1e-6, ... 
    'TolFun',1e-6,... 
    'TolGradCon',1e-6, ...
    'TolX',1e-15); 

% If just 'defaults' passed in, return the default options in X
if isequal(nargin,1) && nargout <= 1 && strcmpi(FUN,'defaults')
    X = defaultopt;
    return
end

if nargin < 11 
    ktrOptsFile = '';
    if nargin < 10 
        options = [];
        if nargin < 9 
            NONLCON = [];
            if nargin < 8 
                UB = [];
                if nargin < 7 
                    LB = [];
                    if nargin < 6 
                        Beq = [];
                        if nargin < 5 
                            Aeq =[];
                            if nargin < 4
                                B = [];
                                if nargin < 3
                                    A = [];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if nargin < 2
    error(message('optim:ktrlink:AtLeastTwoInputs'))
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(X,A,B,Aeq,Beq,LB,UB);
catch ME
    if strcmp(ME.identifier,'MATLAB:datatypes:superiorfloat')
       dataType = 'notDouble';
    end
end

if ~strcmp(dataType,'double')
    error(message('optim:ktrlink:NonDoubleInput'))
end

if nargout > 4
    computeLambda = 1;
else
    computeLambda = 0;
end

% Get maximum integer that will allow us to use linear indexing (platform dependent)
if ~any(strcmp(computer,{'PCWIN','GLNX86'}))
    maxIntIdx = intmax('int64');
else
    maxIntIdx = intmax('int32');
end
    

XOUT = full(X(:));
numberOfVariables = length(XOUT);
% Check for empty X
if isequal(numberOfVariables,0)
    error(message('optim:ktrlink:EmptyX'));
end

% Set to column vectors
B = B(:);
Beq = Beq(:);

if isempty(NONLCON)
    constflag = false;
else
    constflag = true;
end

% Check OPTIMSET created options or KNITRO options file
[knitroOptions,options,hessFcn,hessMultFcn] = ...
    getOpts(options,defaultopt,constflag,ktrOptsFile);

% Make sure bound vectors are the correct length, and bound values to 
% +/-Inf where no bounds exist
[XOUT,lbound,ubound] = checkbounds(XOUT,LB,UB,numberOfVariables);

ktrInf = 1e20;  % KNITRO Inf
% Set Inf values to KNITRO Inf
lbound(isinf(lbound)) = -ktrInf;
ubound(isinf(ubound)) = ktrInf;

gradflag = isequal(knitroOptions.gradopt,1);  % Set gradflag based on KNITRO option gradopt
gradconstflag = gradflag;       % gradconstflag also represented by KNITRO option gradopt

% Convert to function-handle as needed
if ~isempty(FUN)  % Detect empty string, empty matrix, empty cell array
    % lenVarIn and funValCheck hard-coded to 0 and false, respectively, 
    % because neither are supported for objectives
    funfcn = optimfcnchk(FUN,'fmincon',0,false,gradflag);
    objFun = funfcn{3}; % Objective function handle
else
    error(message('optim:ktrlink:InvalidFUN'));
end

if constflag % If there are nonlinear constraints
    % lenVarIn, funValCheck, and hessflag hard-coded to 0, false, and
    % false, respectively, because none are supported for constraints
    confcn = optimfcnchk(NONLCON,'fmincon',0,false,gradconstflag,false,true);
    conFun = confcn{3}; % Nonlinear constraint function handle
end

X(:) = XOUT;

% Evaluate function
GRAD = zeros(numberOfVariables,1);
try
    f = feval(objFun,X);
catch userFcn_ME
    optim_ME = MException('optim:ktrlink:ObjectiveError', ...
        getString(message('optim:ktrlink:ObjectiveError')));
    userFcn_ME = addCause(userFcn_ME,optim_ME);
    rethrow(userFcn_ME)
end

% Check that the objective value is a scalar
if numel(f) ~= 1
    error(message('optim:ktrlink:NonScalarObj'))
end

% Evaluate constraints and initialize gradients to zeros, if necessary
cineq = [];     % Needed to perform size check
ceq = [];
if constflag
    try
        [ctmp,ceqtmp] = feval(conFun,X);
    catch userFcn_ME
        optim_ME = MException('optim:ktrlink:NonlconError', ...
            getString(message('optim:ktrlink:NonlconError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
    cineq = ctmp(:); 
    ceq = ceqtmp(:);
end % if constflag
% Accounting for number of constraints and their derivatives
non_eq = length(ceq);
non_ineq = length(cineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);

cGRAD = sparse(numberOfVariables,non_ineq);  % These get used later in formKtrJacobian
ceqGRAD = sparse(numberOfVariables,non_eq);

tot_lin = lin_eq + lin_ineq;    % Total number of linear constraints
tot_nlin = non_eq + non_ineq;   % Total number of nonlinear constraints
tot_con = tot_lin + tot_nlin;   % Total number of constraints

% Indices and sizes are represented by signed 32-bit integers in KNITRO.
% As a safety check, make sure that these numbers do not overflow here.
if tot_con > intmax('int32') || numberOfVariables > intmax('int32')
    error(message('optim:ktrlink:ProblemTooLarge'));
end
    
% Assign the types of the constraints:
ctype(1:tot_lin) = 1;                       % Linear
ctype(tot_lin+1:(tot_lin+tot_nlin)) = 0;    % Nonlinear

if ~isempty(Aeq) && Aeqcol ~= numberOfVariables
    error(message('optim:ktrlink:WrongNumberOfColumnsInAeq'))
end
if ~isempty(A) && Acol ~= numberOfVariables
    error(message('optim:ktrlink:WrongNumberOfColumnsInA'))
end

% Make sure that linear constraints are consistent A,B,Aeq,Beq
% Row consistency check.
if ~isempty(Aeq) && ~isequal(size(Aeq,1),length(Beq))
        error(message('optim:ktrlink:AeqAndBeqInconsistent'))
end
if ~isempty(A) && ~isequal(size(A,1),length(B))
        error(message('optim:ktrlink:AinAndBinInconsistent'))
end

% Find the non-zeros in the Hessian sparsity pattern and format the data for KNITRO
if isempty(options.HessPattern)
    if isequal(knitroOptions.hessopt,1) % User-provides exact, dense Hessian        
        % Number of elements in upper triangle
        nnzhessPattern = (numberOfVariables*(numberOfVariables + 1))/2; 
        Hrow = zeros(nnzhessPattern,1);
        Hcol = zeros(nnzhessPattern,1);
        idx = 1;
        % Create row and column indices of dense upper-triangular matrix
        % without actually forming the matrix
        for k = 1:numberOfVariables
            Hrow(idx:idx+k-1) = 1:k;
            Hcol(idx:idx+k-1) = k*ones(k,1);
            idx = idx + k;
        end
    else
        Hrow = [];Hcol = [];nnzhessPattern = 0; 
    end
else
    options.HessPattern = triu(options.HessPattern);    % KNITRO only interested in upper triangle
    nnzhessPattern = nnz(options.HessPattern);
    [mHpattern,nHpattern] = size(options.HessPattern);  % Check size of the "pattern" matrix
    if (mHpattern ~= numberOfVariables) || (nHpattern ~= numberOfVariables)
        error(message('optim:ktrlink:InvalidHessPatternDims', numberOfVariables, numberOfVariables));
    end
    [Hrow,Hcol] = find(options.HessPattern);   % Find non-zero locations in upper triangle
end

hessVectorIn = zeros(numberOfVariables,1);  % Create input for solveKTR

% Begin computing constraint Jacobians
maxfunind = 0;    % Index to stack different type of constraint Jacobians
numNonlinConstr = non_ineq + non_eq;
% Flag indicating presence of Jacobian patterns
haveJpattern = ~isempty(options.JacobPattern);

% Initialize indvar,indfun,cjac
nnz_lin = nnz(A);
nnz_lineq = nnz(Aeq);
if constflag
    if haveJpattern
        nnz_non = nnz(options.JacobPattern); % Non-zeros in J pattern matrix
    else
        nnz_non = numNonlinConstr * numberOfVariables; % Total elements in matrix
    end
else
    nnz_non = 0;
end
nnzj = nnz_lin + nnz_lineq + nnz_non;   % Total non-zeros in Jacobian
% Pre-allocate arrays to hold non-zeros and their indices
indvar = zeros(nnzj,1);
indfun = zeros(nnzj,1);
cjac   = zeros(nnzj,1);

% Populate the constraint Jacobian corresponding to linear constraints
if nnz_lin ~= 0
    [indfun,indvar,cjac] = find(A);
    maxfunind = lin_ineq;
end

nnzj = nnz_lin; % Re-assign nnzj to be pointer to current location in cjac
if nnz_lineq ~= 0
    [indfun(nnzj+1:nnzj+nnz_lineq),indvar(nnzj+1:nnzj+nnz_lineq),...
        cjac(nnzj+1:nnzj+nnz_lineq)] = find(Aeq);
    indfun(nnzj+1:nnzj+nnz_lineq) = indfun(nnzj+1:nnzj+nnz_lineq) + maxfunind;
end

% Now the constraint Jacobian corresponding to nonlinear constraints
maxfunind = tot_lin;        % Update row pointer in "virtual" Jacobian array
nnzj = nnz_lin + nnz_lineq; % Re-assign nnzj to be pointer to current location in cjac

if constflag                % If there are nonlinear constraints
    % Get indices of nonzeros for the nonlinear inequalities
    [indfun,indvar,nnzj,Jidx,subIndexJac,Jrow,Jcol] = getJacIdx(haveJpattern, ...
        options.JacobPattern,numNonlinConstr,numberOfVariables,indfun,indvar, ...
        nnzj,nnz_non,maxfunind,maxIntIdx);
    
    % Format constraint Jacobian for KNITRO
    cjac = formKtrJacobian(cGRAD,ceqGRAD,Jidx,haveJpattern,subIndexJac,...
        cjac,nnz_lineq,nnz_lin,nnz_non,Jrow,Jcol);
end

% Convert to C array indexing
indfun = indfun-1;
indvar = indvar-1;

% Initial value of constraint functions
if lin_eq > 0
    c_lin_eq = Aeq*XOUT; 
else
    c_lin_eq = [];
end
if lin_ineq > 0
    c_lin_ineq = A*XOUT; 
else
    c_lin_ineq = [];
end
c = [c_lin_ineq; c_lin_eq; cineq; ceq];

try
    % Initialize KNITRO and get problem identifier and status
    [kcontext,status,knitroOptions.algorithm,Hrow,Hcol] = ...
        setupKTR(numberOfVariables,tot_con,XOUT,lbound,ubound,B,Beq,ctype,...
        lin_eq,lin_ineq,non_eq,nnzj,indvar,indfun,nnzhessPattern,Hrow,Hcol,knitroOptions);
    % Delete problem instance, free memory, etc, when exiting this function.
    onCleanup(@()closeKTR(kcontext));
catch KtrSetupMExcept
    if strcmpi(KtrSetupMExcept.identifier,'MATLAB:invalidMEXFile')
        error(message('optim:ktrlink:MissingLib'));
    else
       rethrow(KtrSetupMExcept); 
    end
end

% Convert row-column indices to linear indices.
hessIdx = sub2ind([numberOfVariables numberOfVariables],Hrow,Hcol);
nnzhessPattern = length(Hrow);
% Protect against integer overflow for large problems
if max(hessIdx) > maxIntIdx
    subIndexHess = true;
    hess = zeros(nnzhessPattern,1);
else
    subIndexHess = false;
    hess = [];
end

% KNITRO specific status flags:
KTR_RC_EVAL_FC = 1; % Evaluate the objective function and constraints
KTR_RC_EVAL_GA = 2; % Evaluate the gradients
KTR_RC_EVAL_H  = 3; % Evaluate the Hessian
KTR_RC_EVAL_HV = 7; % Evaluate the Hessian-vector product

% Initialize loop variables
KnitroDone = false;     % Flag which indicates whether KNITRO is done
checkedFunEval = false; % Indicates 1st gradient eval, needed to ensure we only error check once
checkedGradEval = false;% Indicates 1st gradient eval, needed to ensure we only error check once
checkedHessEval = false;% Indicates 1st HessFcn or HessMult eval, needed to ensure we only error check once
Xinit = XOUT;           % Initial X0 from user, KNITRO may shift XOUT to make feasible wrt bounds
FVAL = f;               % FVAL is output from solveKTR, f is the input equivalent

% Initialize LAMBDA structure, if user provides Hessian
if isequal(knitroOptions.hessopt,1) || isequal(knitroOptions.hessopt,5) || computeLambda
    LAMBDA.ineqlin    = zeros(lin_ineq,1);
    LAMBDA.eqlin      = zeros(lin_eq,1);
    LAMBDA.lower      = zeros(numberOfVariables,1);
    LAMBDA.upper      = zeros(numberOfVariables,1);
    LAMBDA.ineqnonlin = zeros(non_ineq,1);
    LAMBDA.eqnonlin   = zeros(non_eq,1);
end

% Main loop
while ~KnitroDone
    if isequal(status,KTR_RC_EVAL_FC)
        % Evaluate function and constraint values only if KNITRO shifted X0
        % and we haven't evaluated at this point yet
        if checkedFunEval
            f = feval(objFun,X);          % Evaluate objective
            if constflag
                [cineq(:),ceq(:)] = feval(conFun,X);
            end
            % Compute constraints at current X
            if (lin_eq > 0)
                c_lin_eq = Aeq*XOUT;
            end
            if (lin_ineq > 0)
                c_lin_ineq = A*XOUT;
            end
            c = [c_lin_ineq; c_lin_eq; cineq; ceq];
        else
            if any(~isequal(XOUT,Xinit)) % KNITRO shifted X0 to be within bounds, re-evaluate at new X0
                f = feval(objFun,X);
                checkedFunEval = true;
                if constflag
                    [cineq(:),ceq(:)] = feval(conFun,X);
                end

                if lin_eq > 0
                    c_lin_eq = Aeq*XOUT;
                end
                if lin_ineq > 0
                    c_lin_ineq = A*XOUT;
                end
                c = [c_lin_ineq; c_lin_eq; cineq; ceq];
            end
        end
    end

    if isequal(status,KTR_RC_EVAL_GA)     % Evaluate user-supplied gradients 
        if checkedGradEval
            [~,GRAD(:)] = feval(objFun,X);

            if constflag
                [cineq(:),ceq(:),cGRAD,ceqGRAD] = feval(conFun,X);

                cjac = formKtrJacobian(cGRAD,ceqGRAD,Jidx,haveJpattern,subIndexJac,...
                                       cjac,nnz_lineq,nnz_lin,nnz_non,Jrow,Jcol);
            end
        else                        % Perform some error checking
            try
                [~,GRAD(:)] = feval(objFun,X);
            catch userFcn_ME
                optim_ME = MException('optim:ktrlink:ObjGradError', ...
                    getString(message('optim:ktrlink:ObjGradError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            if constflag            % If we have nonlinear constraints
                try
                    [cineq(:),ceq(:),cGRAD,ceqGRAD] = feval(conFun,X);
                catch userFcn_ME
                    optim_ME = MException('optim:ktrlink:ConstrGradError', ...
                        getString(message('optim:ktrlink:ConstrGradError')));
                    userFcn_ME = addCause(userFcn_ME,optim_ME);
                    rethrow(userFcn_ME)
                end
                % Check for consistency in constraint derivatives
                [cgrow, cgcol] = size(cGRAD);
                [ceqgrow, ceqgcol] = size(ceqGRAD);
                if  (~isequal(non_ineq,0) && cgrow ~= numberOfVariables) || cgcol ~= non_ineq
                    error(message('optim:ktrlink:WrongSizeGradNonlinIneq'))
                end
                if (~isequal(non_eq,0) && ceqgrow ~= numberOfVariables) || ceqgcol ~= non_eq
                    error(message('optim:ktrlink:WrongSizeGradNonlinEq'))
                end
                % Format Jacobian for KNITRO
                cjac = formKtrJacobian(cGRAD,ceqGRAD,Jidx,haveJpattern,subIndexJac,...
                    cjac,nnz_lineq,nnz_lin,nnz_non,Jrow,Jcol);
            end
            checkedGradEval = true;
        end
        GRAD = full(GRAD);  % Enforce that GRAD is a dense vector as needed by KNITRO
    end % Gradient calculation request

    if isequal(status,KTR_RC_EVAL_H)
        % Read Lagrange multipliers from KNITRO, assuming that KNITRO will
        % populate the LMBDOUT vector.
        LAMBDA.ineqnonlin = LMBDOUT(tot_lin+1:tot_lin+non_ineq);   
        LAMBDA.eqnonlin   = LMBDOUT(tot_lin+non_ineq+1:tot_con);
        if checkedHessEval
            HESS = feval(hessFcn,X,LAMBDA);
        else            
            try % Evaluating Hessian Function
                HESS = feval(hessFcn,X,LAMBDA);
            catch userFcn_ME
                optim_ME = MException('optim:ktrlink:HessFcnError', ...
                    getString(message('optim:ktrlink:HessFcnError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            checkedHessEval = true;
        end
        nnzHESS = nnz(triu(HESS));
        if nnzHESS > nnzhessPattern
            error(message('optim:ktrlink:InvalidHessStructure'));
        end
        % HESS is the Hessian as computed by the user, in matrix form.
        % hess is only the non-zeros of the Hessian requested by KNITRO in vector format.
        % Extract elements from user Hessian according to pattern
        if ~subIndexHess        
            hess = full(HESS(hessIdx));   
        else
            for k = 1:nnzhessPattern
                hess(k) = HESS(Hrow(k),Hcol(k));
            end
        end
    end
    
    if isequal(status,KTR_RC_EVAL_HV)
        % Read Lagrange multipliers from KNITRO
        LAMBDA.ineqnonlin = LMBDOUT(tot_lin+1:tot_lin+non_ineq);
        LAMBDA.eqnonlin   = LMBDOUT(tot_lin+non_ineq+1:tot_con);
        if checkedHessEval
            % hessVectorIn is the result of the Hessian-vector product.
            % hessVectorOut is the vector with which the Hessian is to be multiplied.
            hessVectorIn = feval(hessMultFcn,X,LAMBDA,hessVectorOut);
        else
            try % Evaluating Hessian-vector product
                hessVectorIn = feval(hessMultFcn,X,LAMBDA,hessVectorOut);
            catch userFcn_ME
                optim_ME = MException('optim:ktrlink:HessMultError', ...
                    getString(message('optim:ktrlink:HessMultError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            checkedHessEval = true;
        end
    end        
    % Flag indicating that the evaluation was successful
    evalStatus = 0;

    % Call KNITRO solver
    [status,XOUT,FVAL,LMBDOUT,hessVectorOut] = ...
        solveKTR(kcontext,numberOfVariables,tot_con,XOUT,evalStatus,f,c,GRAD,cjac,hess,hessVectorIn);
    % Assign X and f to the output from KNITRO
    X(:) = XOUT;
    f = FVAL;
    KnitroDone = status <= 0;
end % Main loop

EXITFLAG = status;  % status flag also serves as KNITRO exitflag

% Read firstorderopt, iterations, and maxconstr from KNITRO
OUTPUT.firstorderopt   = getSolKTR(kcontext,'abs_opt_error');
OUTPUT.iterations      = getSolKTR(kcontext,'major_iters');
OUTPUT.constrviolation = getSolKTR(kcontext,'abs_feas_error');
OUTPUT.funcCount       = getSolKTR(kcontext,'number_FC_evals');

% Updated with information from KNITRO
switch knitroOptions.algorithm
    case 0
        OUTPUT.algorithm = 'KNITRO Auto-select';
    case 1
        OUTPUT.algorithm = 'KNITRO Interior/Direct';
    case 2
        OUTPUT.algorithm = 'KNITRO Interior/CG';
    case 3
        OUTPUT.algorithm = 'KNITRO Active-Set';
end

% Fill in MATLAB structure of Lagrange multipliers, LAMBDA, with values
% returned from KNITRO.
if computeLambda
    LAMBDA.ineqlin    = LMBDOUT(1:lin_ineq);
    LAMBDA.eqlin      = LMBDOUT(lin_ineq+1:tot_lin);
    LAMBDA.ineqnonlin = LMBDOUT(tot_lin+1:tot_lin+non_ineq);
    LAMBDA.eqnonlin   = LMBDOUT(tot_lin+non_ineq+1:tot_con);
    % Determine which multipliers correspond to upper or lower bounds.
    if isempty(LB)
        LAMBDA.upper = LMBDOUT(tot_con+1:end);
    elseif isempty(UB)
        LAMBDA.lower = LMBDOUT(tot_con+1:end);
    else % Problem has lower and upper bounds
        for i = 1:numberOfVariables
            if abs(X(i) - lbound(i)) < abs(X(i) - ubound(i))
                LAMBDA.lower(i) = LMBDOUT(tot_con+i);
            else
                LAMBDA.upper(i) = LMBDOUT(tot_con+i);
            end
        end
    end
end % if computeLambda

%-------------------------------------------------------------------------------
function [knitroOptions,options,hessFcn,hessMultFcn] = ...
         getOpts(options,defaultopt,constflag,ktrOptsFile)
%getOpts a utility function that extracts and checks all KTRLINK-based options
%   getOpts extracts all fmincon-based options relevant to ktrlink and
%   returns the necessary variables used by KNITRO and the rest of
%   ktrlink.

% Check HessFcn, HessMult, HessPattern, and JacobPattern here, since these
% are needed when there is a KNITRO options file or a full OPTIONS
% structure.
hessFcn = optimget(options,'HessFcn',defaultopt,'fast');
hessMultFcn = optimget(options,'HessMult',defaultopt,'fast');

hessPattern = optimget(options,'HessPattern',defaultopt,'fast'); % Hessian sparsity pattern
if ischar(hessPattern)
    if strcmpi(hessPattern,'sparse(ones(numberofvariables))')
        options.HessPattern = [];
    else
        error(message('optim:ktrlink:InvalidHessPattern'))
    end
end

JacPattern = optimget(options,'JacobPattern',defaultopt,'fast'); % Hessian sparsity pattern
if ischar(JacPattern)
    if strcmpi(JacPattern,'sparse(ones(Jrows,Jcols))')
        options.JacobPattern = [];
    else
        error(message('optim:ktrlink:InvalidJacobPattern'))
    end
end

% Check here to see if user has KNITRO options file
if isempty(ktrOptsFile)  % Proceed checking all options via OPTIMGET
    knitroOptions.ktroptsfile = '';
    switch optimget(options,'Display',defaultopt,'fast')
        case {'off','none'}
            knitroOptions.outlev = 0;
        case {'iter','iter-detailed'}
            knitroOptions.outlev = 4;
        case {'notify','notify-detailed'}
            knitroOptions.outlev = 1;
            warning(message('optim:ktrlink:NotifyUnavailable'));
        case {'final','final-detailed'}
            knitroOptions.outlev = 1;
        otherwise
            knitroOptions.outlev = 1;
    end

    % Find out what algorithm user wants to run, re-map to KNITRO code:
    knitro_algorithm = optimget(options,'Algorithm',defaultopt,'fast');
    subProbMethod = optimget(options,'SubproblemAlgorithm',defaultopt,'fast');
    if iscell(subProbMethod)
        if isequal(numel(subProbMethod),2)          % If 2 element cell-array, pivot tolerance given
            knitroOptions.tolpivot = subProbMethod{2};  % Pivot tolerance for LDL factorization
        end
        subProbMethod = subProbMethod{1};
    else
        knitroOptions.tolpivot = 1e-8;              % KNITRO default for pivot tolerance
    end

    switch knitro_algorithm
        case 'interior-point'
            if strcmp(subProbMethod,'ldl-factorization')  
                knitroOptions.algorithm = 1;  % KNITRO algorithm 1 is interior-point using direct factorization
            elseif strcmp(subProbMethod,'cg')
                knitroOptions.algorithm = 2;  % KNITRO algorithm 2 is interior-point using CG
            else
                error(message('optim:ktrlink:BadSubproblemAlg'));
            end
        case 'active-set'
            knitroOptions.algorithm = 3;
        otherwise
            error(message('optim:ktrlink:Algorithm', knitro_algorithm));
    end
    
    % Check which gradients are given (objective function, constraints, both,
    % none). KNITRO requires both to be computed by the same method. 
    gradflag = strcmp(optimget(options,'GradObj',defaultopt,'fast'),'on');
    gradconstflag = strcmp(optimget(options,'GradConstr',defaultopt,'fast'),'on');
    finDiffMethod = optimget(options,'FinDiffType',defaultopt,'fast');
    % If gradients are given by different methods, make them the same and warn
    if xor(gradflag,gradconstflag)
        if constflag
            gradflag = false; 
            warning(message('optim:ktrlink:GradMethodInconsistent'));
        end
    end
    % gradflag and gradconstflag should be the same now, so only check gradflag.
    if gradflag  % user-provided
        knitroOptions.gradopt = 1;
    else
        switch finDiffMethod
            case 'forward'
                knitroOptions.gradopt = 2;
            case 'central'
                knitroOptions.gradopt = 3;
        end
    end % if gradflag

    knitroOptions.maxit = optimget(options,'MaxIter',defaultopt,'fast');
    knitroOptions.xtol = optimget(options,'TolX',defaultopt,'fast');
    knitroOptions.feastol = optimget(options,'TolCon',defaultopt,'fast');
    knitroOptions.scale = strcmpi(optimget(options,'ScaleProblem',defaultopt,'fast'),'obj-and-constr');
    knitroOptions.opttol = optimget(options,'TolFun',defaultopt,'fast');

    if strcmpi(optimget(options,'AlwaysHonorConstraints',defaultopt,'fast'),'bounds')
        knitroOptions.honorbounds = 1;    % All iterations satisfy bounds
    else
        knitroOptions.honorbounds = 0;    % Do not have to satisfy bounds
    end

    initDelta = optimget(options,'InitTrustRegionRadius',defaultopt,'fast');
    if ischar(initDelta)
        if strcmpi(initDelta,'sqrt(numberofvariables)')
            knitroOptions.delta = 1;    % Set to KNITRO default
        end
    else
        knitroOptions.delta = initDelta;
    end

    knitroOptions.maxcgit = optimget(options,'MaxProjCGIter',defaultopt,'fast'); % Get Proj. CG iterations limit
    if ischar(knitroOptions.maxcgit)
        if strcmpi(knitroOptions.maxcgit,'2*(numberofvariables-numberofequalities)')
            knitroOptions.maxcgit = 0;  % KNITRO adjusts this quantity automatically
        end
    end
    
    knitroOptions.initmu = optimget(options,'InitBarrierParam',defaultopt,'fast');
    knitroOptions.objrange = optimget(options,'ObjectiveLimit',defaultopt,'fast');
    
    hessType = optimget(options,'Hessian',defaultopt,'fast'); % which of 5 types of Hessians
    if iscell(hessType)
        if isequal(numel(hessType),2)       % % If 2 element cell-array, LBFGS pair limit given
            knitroOptions.lmsize = hessType{2};
        end
        hessType = hessType{1};
    else
        knitroOptions.lmsize = 10;       % default KNITRO value
    end

    knitroOptions.hessopt = 2;      % BFGS by default, for 'bfgs','on', and 'off'
    switch hessType
        case 'user-supplied'
            if ~isempty(hessFcn)
                knitroOptions.hessopt = 1;    % User Provided Exact Hessian
            elseif ~isempty(hessMultFcn)
                % Check if name clash
                functionNameClashCheck('HessMult',hessMultFcn, ...
                    'hessMult_optimInternal','optim:ktrlink:HessMultNameClash');
                knitroOptions.hessopt = 5;
            else % user-supplied selected, but no function
                error(message('optim:ktrlink:HessUserFcnNotGiven'));
            end
        case 'lbfgs'
            knitroOptions.hessopt = 6;
        case 'fin-diff-grads'
            knitroOptions.hessopt = 4;
    end
else
    % User has specified a KNITRO-options file. Relevant fields of OPTIONS:
    % HessFcn, HessMult, HessPattern, and JacobPattern already checked
    % above.
    
    % Create local structure corresponding to KNITRO options needed by
    % the rest of ktrlink. In this case, most options will be read in directly by
    % KNITRO from the options text file and therefore we don't need to set them.
    
    % We do not know what Hessian setting the user requests inside the
    % options text file. We will check hessopt inside setupKTR, but must
    % make sure that we use the HessPattern if one has been provided.
    if isempty(options.HessPattern)
        % In this case, set hessopt to empty and read the actual hessopt
        % value after letting KNITRO read the text file. At that point we
        % will create a dense upper triangular HessPattern if hessopt = 1
        % (exact).
        hessopt = [];
    else
        % In this case, set hessopt to 1 (exact) so that we read the user's
        % HessPattern. the pattern and pass it to KNITRO. It doesn't matter
        % if the user is not providing an exact Hessian. KNITRO will just
        % ignore the pattern if it doesn't need it. 
        hessopt = 1;
    end
    
    knitroOptions = struct('gradopt',2,'hessopt',hessopt,'ktroptsfile',ktrOptsFile);
end  % if isempty(ktrOptsFile)

%-------------------------------------------------------------------------------
function [indfun,indvar,nnzj,JacIdx,subIndexJac,Jrow,Jcol] = getJacIdx(havePattern,...
    JacStrct,numConstr,numVars,indfun,indvar,nnzj,nnzCon,maxfunind,maxIntIdx)
%getJacIdx A utility function of ktrlink that finds the constraint Jacobian indices
%   getJacIdx takes the Jacobian sparsity structure and finds the subscript indices
%   of the non-zero elements for KNITRO, as well as returning the linear
%   indices for extracting the non-zeros of the computed Jacobian later.
   
if havePattern
    [mStr,nStr] = size(JacStrct);   % Check size of the pattern matrix
    if (mStr ~= numConstr) || (nStr ~= numVars)
        error(message('optim:ktrlink:InvalidJacStrSize', sprintf( '%g', numConstr ), sprintf( '%g', numVars )));
    end

    [indfun(nnzj+1:nnzj+nnzCon),indvar(nnzj+1:nnzj+nnzCon)] = find(JacStrct);
    % Convert to linear index
    JacIdx = sub2ind([mStr,nStr],indfun(nnzj+1:nnzj+nnzCon),indvar(nnzj+1:nnzj+nnzCon));
    
    % Protect against integer overflow for large problems
    if max(JacIdx) > maxIntIdx
        subIndexJac = true;
        Jrow = indfun(nnzj+1:nnzj+nnzCon);
        Jcol = indvar(nnzj+1:nnzj+nnzCon);
    else
        subIndexJac = false;
        Jrow = [];Jcol = [];
    end
    
    % Offset row indices by number of existing rows in "virtual" Jacobian array
    indfun(nnzj+1:nnzj+nnzCon) = indfun(nnzj+1:nnzj+nnzCon) + maxfunind;
else
    indfun(nnzj+1:nnzj+nnzCon) = repmat((1:numConstr)',numVars,1) + maxfunind;
    tmpMatrix = repmat((1:numVars),numConstr,1);
    indvar(nnzj+1:nnzj+nnzCon) = tmpMatrix(:);
    % Needed as input for formKtrJacobian
    JacIdx = [];Jrow = [];Jcol = [];
    subIndexJac = false;
end
nnzj = nnzj + nnzCon;                % Add to counter of non-zeros

%-------------------------------------------------------------------------------
function cjac = formKtrJacobian(cGRAD,ceqGRAD,Jidx,haveJpattern,subIndexJac,...
                                cjac,nnz_lineq,nnz_lin,nnz_non,Jrow,Jcol)
%formKtrJacobian A utility function that forms the constraint Jacobian for KNITRO
%   formKtrJacobian takes the computed constraint gradients and formats
%   them for use by KNITRO. The formatted Jacobian is returned in vector
%   cjac. In addition, formKtrJacobian returns the count of non-zeros in
%   the constraint Jacobian.

% Set counter of non-zeros back to the beginning of the nonlinear constraints
nnzj = nnz_lineq + nnz_lin; 
JAC = [cGRAD'; ceqGRAD']; % Form Jacobian to match shape of sparsity structure
if haveJpattern
    nnzJac = nnz(JAC);       % Non-zeros in computed Jacobian
    if nnzJac > nnz_non
        error(message('optim:ktrlink:formKtrJacobian:InvalidJacStructure'));
    end
    % extract elements from user Jacobian according to pattern
    if ~subIndexJac
        cjac(nnzj+1:nnzj+nnz_non) = full(JAC(Jidx));
    else
        for k = 1:nnz_non
           cjac(nnzj+k) = JAC(Jrow(k),Jcol(k));
        end
    end
else
    cjac(nnzj+1:nnzj+nnz_non) = JAC;
end

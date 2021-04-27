function [x,FVAL,EXITFLAG,OUTPUT,LAMBDA] = fseminf(FUN,x,ntheta,SEMINFCON,A,B,Aeq,Beq,LB,UB,options,varargin) 
%FSEMINF solves semi-infinite constrained optimization problems.
%   FSEMINF attempts to solve problems of the form:
%
%          min { F(x) | C(x) <= 0 , Ceq(x) = 0 , PHI(x,w) <= 0 }
%           x
%   for all w in an interval. 
%
%   X = FSEMINF(FUN,X0,NTHETA,SEMINFCON) starts at X0 and finds minimum X 
%   to the function FUN constrained by NTHETA semi-infinite constraints in 
%   the function SEMINFCON. FUN accepts vector input X and returns the 
%   scalar function value F evaluated at X. Function SEMINFCON accepts 
%   vector inputs X and S and return a vector C of nonlinear inequality
%   constraints, a vector Ceq of nonlinear equality constraints and
%   NTHETA semi-infinite inequality constraint matrices, PHI_1, PHI_2, ...,
%   PHI_NTHETA, evaluated over an interval. S is a recommended sampling
%   interval, which may or may not be used. 
%
%   X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B) also tries to satisfy the  
%   linear inequalities A*X <= B.
%
%   X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq) minimizes subject to 
%   the linear equalities Aeq*X = Beq as well.  (Set A=[] and B=[] if no
%   inequalities exist.)
%
%   X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB) defines a set of
%   lower and upper bounds on the design variables, X, so that the
%   solution is in the range LB <= X <= UB. Use empty matrices for LB and U
%   if no bounds exist.  Set LB(i) = -Inf if X(i) is unbounded below; set
%   UB(i) = Inf if X(i) is unbounded above.
%
%   X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB,OPTIONS) 
%   minimizes with the default optimization parameters replaced by values 
%   in the structure OPTIONS, an argument created with the OPTIMSET 
%   function. See OPTIMSET for details. Used options are Display, TolX, 
%   TolFun, TolCon, DerivativeCheck, Diagnostics, FunValCheck, GradObj, 
%   MaxFunEvals, MaxIter, DiffMinChange, DiffMaxChange, PlotFcns,  
%   OutputFcn, and TypicalX. Use the GradObj option to specify that FUN may
%   be called with two output arguments where the second, G, is the partial 
%   derivatives of the function df/dX, at the point X:  
%   [F,G] = feval(FUN,X). 
%
%   X = FSEMINF(PROBLEM) solves the semi-infinite constrained problem 
%   defined in PROBLEM. PROBLEM is a structure with the function FUN in 
%   PROBLEM.objective, the start point in PROBLEM.x0, the number of 
%   semi-infinite constraints in PROBLEM.ntheta, the nonlinear and 
%   semi-infinite constraint function in PROBLEM.seminfcon, the linear 
%   inequality constraints in PROBLEM.Aineq and PROBLEM.bineq, the linear 
%   equality constraints in PROBLEM.Aeq and PROBLEM.beq, the lower bounds 
%   in PROBLEM.lb, the upper bounds in PROBLEM.ub, the options structure 
%   in PROBLEM.options, and solver name 'fseminf' in PROBLEM.solver. Use 
%   this syntax to solve at the command line a problem exported from 
%   OPTIMTOOL. The structure PROBLEM must have all the fields. 
%
%   [X,FVAL] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) returns the value of 
%   the objective function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) returns an 
%   EXITFLAG that describes the exit condition of FSEMINF. Possible values  
%   of EXITFLAG and the corresponding exit conditions are listed below. See
%   the documentation for a complete description.
%
%     1  FSEMINF converged to a solution.
%     4  Computed search direction too small.
%     5  Predicted change in objective function too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  No feasible point found.
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) returns 
%   a structure OUTPUT with the number of iterations taken in  
%   OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the norm of the final step in OUTPUT.stepsize, the 
%   final line search steplength in OUTPUT.lssteplength, the algorithm used
%   in OUTPUT.algorithm, the first-order optimality in  
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...)
%   returns the Lagrange multiplier at the solution X: LAMBDA.lower for
%   LB, LAMBDA.upper for UB, LAMBDA.ineqlin for the linear inequalities,
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear
%   equalities.
%
%   Examples
%     FUN and SEMINFCON can be specified using @:
%        x = fseminf(@myfun,[2 3 4],3,@myseminfcon)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = x(1)*cos(x(2))+ x(3)^3;
%
%   and myseminfcon is a MATLAB function such as:
%
%       function [C,Ceq,PHI1,PHI2,PHI3,S] = myseminfcon(X,S)
%       C = [];     % Code to compute C and Ceq: could be 
%                   % empty matrices if no constraints.
%       Ceq = [];
%       if isnan(S(1,1))
%          S = [...] ; % S has ntheta rows and 2 columns
%       end
%       PHI1 = ... ;       % code to compute PHI's
%       PHI2 = ... ;
%       PHI3 = ... ;
%
%   See also OPTIMSET, @, FGOALATTAIN, LSQNONLIN.

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.22 $  $Date: 2011/06/30 16:36:46 $
%       

defaultopt = struct( ...
    'DerivativeCheck','off',...
    'Diagnostics','off',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0, ...
    'Display','final',...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward', ...
    'FunValCheck','off',...
    'GradConstr',[],...
    'GradObj','off',...
    'Hessian',[],...  % Hessian and GradConstr not used
    'MaxFunEvals','100*numberOfVariables',...
    'MaxIter',400,...
    'MaxSQPIter','10*max(numberOfVariables,numberOfInequalities+numberOfBounds)',...
    'NoStopIfFlatInfeas','off',...
    'OutputFcn',[],...
    'PhaseOneTotalScaling','off',...
    'PlotFcns',[], ...
    'RelLineSrchBnd',[],...
    'RelLineSrchBndDuration',1,...
    'TolCon',1e-6,...
    'TolConSQP',1e-6,...
    'TolFun',1e-4,...
    'TolX',1e-4,...
    'TypicalX','ones(numberOfVariables,1)'...
    );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 11
    options = [];
    if nargin < 10
        UB = [];
        if nargin < 9
            LB = [];
            if nargin < 8
                Beq = [];
                if nargin < 7
                    Aeq = [];
                    if nargin < 6
                        B = [];
                        if nargin < 5
                            A = [];
                        end
                    end
                end
            end
        end
    end
end

% Detect problem structure input
problemInput = false;
if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,x,ntheta,SEMINFCON,A,B,Aeq,Beq,LB,UB,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fseminf:InputArg'));
    end
end

if nargin < 4 && ~problemInput
    error(message('optim:fseminf:NotEnoughInputs'))
end,

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(x,ntheta,A,B,Aeq,Beq,LB,UB);
catch ME
    if strcmp(ME.identifier,'MATLAB:datatypes:superiorfloat')
       dataType = 'notDouble';
    end
end

if ~strcmp(dataType,'double')
    error(message('optim:fseminf:NonDoubleInput'))
end

initVals.xOrigShape = x;
[sizes.xRows,sizes.xCols] = size(x);
xnew = x(:);
sizes.nVar = length(xnew);

diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');

display = optimget(options,'Display',defaultopt,'fast');
flags.detailedExitMsg = ~isempty(strfind(display,'detailed'));
switch display
    case {'off','none'}
        verbosity = 0;
    case {'notify','notify-detailed'}
        verbosity = 1;
    case {'final','final-detailed'}
        verbosity = 2;
    case {'iter','iter-detailed'}
        verbosity = 3;
    otherwise
        verbosity = 2;
end

% Set to column vectors
B = B(:);
Beq = Beq(:);

[xnew,l,u,msg] = checkbounds(xnew,LB,UB,sizes.nVar);
if ~isempty(msg)
    EXITFLAG = -2;
    [FVAL,LAMBDA] = deal([]);
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.stepsize = [];
    OUTPUT.lssteplength = [];
    OUTPUT.algorithm = 'semi-infinite, SQP, Quasi-Newton, line_search';
    OUTPUT.firstorderopt = [];
    OUTPUT.constrviolation = [];
    OUTPUT.message = msg;
    x(:) = xnew(1:sizes.nVar);
    if verbosity > 0
        disp(msg)
    end
    return
end

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
DerivativeCheck = strcmpi(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');

% Check options needed for finitedifferences
% Write checked DiffMaxChange, DiffMinChage, FinDiffType, FinDiffRelStep,
% GradObj and GradConstr options back into struct for later use
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
if options.DiffMinChange >= options.DiffMaxChange
    error(message('optim:fseminf:DiffChangesInconsistent', sprintf( '%0.5g', options.DiffMinChange ), sprintf( '%0.5g', options.DiffMaxChange )))
end
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx;
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast');
options = validateFinDiffRelStep(sizes.nVar,options,defaultopt);

options.GradObj = optimget(options,'GradObj',defaultopt,'fast');

flags.grad = strcmpi(options.GradObj,'on');
flags.hess = strcmp(optimget(options,'Hessian',defaultopt,'fast'),'on');
if flags.hess
   warning(message('optim:fseminf:UserHessNotUsed'))
   flags.hess = 0;
end

if isempty(SEMINFCON)
   userconstflag = 0;
else
   userconstflag = 1;
end
flags.gradconst = strcmp(optimget(options,'GradConstr',defaultopt,'fast'),'on');
options.GradConstr = 'off';
if flags.gradconst
   warning(message('optim:fseminf:UserConstrGradNotUsed'))
   flags.gradconst = 0;
end

flags.meritFunction = 5;  % formerly options(7)
lenVarIn = length(varargin);
line_search = 1;
% semicon also needs ntheta and FunStr
semargs = 2;

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array   
   funfcn = optimfcnchk(FUN,'fseminf',lenVarIn,funValCheck,flags.grad,flags.hess);
else
   error(message('optim:fseminf:InvalidFUN'))
end

if ~isempty(SEMINFCON)
    % SEMINFCON cannot be an inline since variable number of output arguments
    % Use optimfcnchk to figure out if SEMINFCON is an inline or expression
    [userconfcn, msg] = optimfcnchk(SEMINFCON,'fseminf',lenVarIn,funValCheck,false,false,false,[],ntheta);
    if isa(SEMINFCON,'inline') % an inline object
        error(message('optim:fseminf:NoInlineSEMINFCON'))
    elseif isa(userconfcn,'inline') % an expression turned into an inline by fcnchk
        error(message('optim:fseminf:NoExprSEMINFCON'))
    elseif ~isempty(msg)
        error(message('optim:fseminf:InvalidSEMINFCON'))
    end
   %  semicon is actually called directly by nlconst, but to be
   %   compatible with the other calls to nlconst, call optimfcnchk.
   % Pass 'false' for funValCheck argument as we don't need to check this call.
   [confcn, msg] = optimfcnchk(@semicon,'fseminf',lenVarIn+semargs,false,flags.gradconst,false,true);
   if ~isempty(msg)
      error(message(msg.identifier))
   end
else
   error(message('optim:fseminf:NoSeminfConstr'))
end

if ntheta < 1
   error(message('optim:fseminf:InvalidNumOfConstr'))
end

lenvlb = length(l);
lenvub = length(u);

i = 1:lenvlb;
lindex = xnew(i) < l(i);
if any(lindex)
    xnew(lindex) = l(lindex) + 1e-4;
end
i = 1:lenvub;
uindex = xnew(i) > u(i);
if any(uindex)
    xnew(uindex) = u(uindex);
end
x(:) = xnew;

if xnew(i) < l(i)
    xnew(i) = l(i) + 1e-4;
end
if xnew(i) > u(i)
    xnew(i) = u(i);
end

x(:) = xnew; s = NaN; initVals.g = zeros(sizes.nVar,1); initVals.H = [];
% Evaluate user function to get number of function values at x
switch funfcn{1}
    case 'fun'
        try
            initVals.f = feval(funfcn{3},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:ObjectiveError', ...
                getString(message('optim:fseminf:ObjectiveError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
    case 'fungrad'
        try
            [initVals.f,initVals.g(:)] = feval(funfcn{3},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:ObjectiveError', ...
                getString(message('optim:fseminf:ObjectiveError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
    case 'fungradhess'
        try
            [initVals.f,initVals.g(:)] = feval(funfcn{3},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:ObjectiveError', ...
                getString(message('optim:fseminf:ObjectiveError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
    case 'fun_then_grad'
        try
            initVals.f = feval(funfcn{3},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:ObjectiveError', ...
                getString(message('optim:fseminf:ObjectiveError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        try
            initVals.g(:) = feval(funfcn{4},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:GradientError', ...
                getString(message('optim:fseminf:GradientError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
    case 'fun_then_grad_then_hess'
        try
            initVals.f = feval(funfcn{3},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:ObjectiveError', ...
                getString(message('optim:fseminf:ObjectiveError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
        try
            initVals.g(:) = feval(funfcn{4},x,varargin{:});
        catch userFcn_ME
            optim_ME = MException('optim:fseminf:GradientError', ...
                getString(message('optim:fseminf:GradientError')));
            userFcn_ME = addCause(userFcn_ME,optim_ME);
            rethrow(userFcn_ME)
        end
    otherwise
        error(message('optim:fseminf:InvalidCalltype'))
end
initVals.f = initVals.f(:);
K = cell(1,ntheta);

try
    switch userconfcn{1}
    case 'fun'
        % Though the last two output arguments are not used, here we want to 
        % exercise the user function with its complete list of arguments.
        [ctmp,ceqtmp,K{:},stmp] = feval(userconfcn{3},x,s,varargin{:});
    otherwise
        error(message('optim:fseminf:InvalidCalltype'))
    end
catch userFcn_ME
    optim_ME = MException('optim:fseminf:BadUserFcn', ...
        getString(message('optim:fseminf:BadUserFcn')));
    userFcn_ME = addCause(userFcn_ME,optim_ME);
    rethrow(userFcn_ME)
end
initVals.ncineq = ctmp(:);
initVals.nceq = ceqtmp(:);

sizes.mNonlinEq = length(initVals.nceq);
initVals.gnceq = zeros(sizes.nVar,sizes.mNonlinEq);

% just_user_constraints contains the number of user-provided finite
% nonlinear inequality constraints. The nonlinear inequalities
% corresponding to the semi-infinite constraints will be appended
% to initVals.ncineq after the next call to semicon.
just_user_constraints = length(initVals.ncineq);

POINT =[];
NEWLAMBDA =[];
LAMBDA = [];
FLAG = 2;
OLDLAMBDA = [];
startnlineq = []; % Start index in LAMBDA for nonlinear inequality constraints
switch confcn{1}
    case 'fun'
        % Call semicon to populate the nonlinear inequalities including the
        % ones corresponding to the semi-inifinite constraints.
        ctmp = feval(confcn{3},xnew,LAMBDA,NEWLAMBDA,OLDLAMBDA,...
            POINT,FLAG,s,startnlineq,ntheta,userconfcn,varargin{:});
        initVals.ncineq = ctmp(:);
        initVals.gnc = zeros(sizes.nVar,length(initVals.ncineq));
    otherwise
        error(message('optim:fseminf:UndefinedCalltype'))
end

sizes.mNonlinIneq = length(initVals.ncineq);

% Make sure empty constraint and their derivatives have correct sizes (not 0-by-0):
if isempty(initVals.ncineq)
    initVals.ncineq = reshape(initVals.ncineq,0,1);
end
if isempty(initVals.nceq)
    initVals.nceq = reshape(initVals.nceq,0,1);
end
if isempty(Aeq)
    Aeq = reshape(Aeq,0,sizes.nVar);
    Beq = reshape(Beq,0,1);
end
if isempty(A)
    A = reshape(A,0,sizes.nVar);
    B = reshape(B,0,1);
end

[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);
[cgrow, cgcol]= size(initVals.gnc);
[ceqgrow, ceqgcol]= size(initVals.gnceq);

if Aeqcol ~= sizes.nVar
    error(message('optim:fseminf:InvalidSizeOfAeq', sizes.nVar))
end
if Acol ~= sizes.nVar
    error(message('optim:fseminf:InvalidSizeOfA', sizes.nVar))
end
if  cgrow ~= sizes.nVar || cgcol ~= sizes.mNonlinIneq
    error(message('optim:fseminf:InvalidSizeOfGC', sizes.nVar, sizes.mNonlinIneq))
end
if ceqgrow ~= sizes.nVar || ceqgcol ~= sizes.mNonlinEq
    error(message('optim:fseminf:InvalidSizeOfGCeq', sizes.nVar, sizes.mNonlinEq))
end

OUTPUT.algorithm = 'semi-infinite, SQP, Quasi-Newton, line_search';  % override nlconst output
if diagnostics
    % Do diagnostics on information so far
    diagnose('fseminf',OUTPUT,flags.grad,flags.hess,userconstflag,flags.gradconst,...
        x,sizes.mNonlinEq,just_user_constraints,lin_eq,lin_ineq,LB,UB,funfcn,confcn);
end

% Create default structure of flags for finitedifferences:
% This structure will (temporarily) ignore some of the features that are
% algorithm-specific (e.g. scaling and fault-tolerance) and can be turned
% on later for the main algorithm.
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward');
finDiffFlags.scaleObjConstr = false; % No scaling for now
finDiffFlags.chkFunEval = false;     % No fault-tolerance yet
finDiffFlags.chkComplexObj = false;  % No need to check for complex values
finDiffFlags.isGrad = true;          % Scalar objective
finDiffFlags.hasLBs = isfinite(l);   % Finite lower bounds
finDiffFlags.hasUBs = isfinite(u);   % Finite upper bounds

% Check objective derivatives if user wants (constraint derivatives are not allowed)
if DerivativeCheck && flags.grad
    validateFirstDerivatives(funfcn,[],x, ...
        l,u,options,finDiffFlags,sizes,varargin{:});
end

initVals.H = [];
problemInfo = []; % No problem related data
[x,FVAL,LAMBDA,EXITFLAG,OUTPUT]=...
   nlconst(funfcn,x,l,u,full(A),B,full(Aeq),Beq,confcn,options,defaultopt, ...
   finDiffFlags,verbosity,flags,initVals,problemInfo,ntheta,userconfcn,varargin{:});

if ~isempty(LAMBDA)
   LAMBDA.ineqnonlin = LAMBDA.ineqnonlin(1:just_user_constraints);
end

% Unimplemented feature: multipliers of semi-inf constraints 
% LAMBDA.semi_infinite = lambda.ineqnonlin(just_user_constraints+1:end);
OUTPUT.algorithm = 'semi-infinite, SQP, Quasi-Newton, line_search';  % override nlconst output

% end seminf

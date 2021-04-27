function [x,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = fgoalattain(FUN,x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FGOALATTAIN solves the multi-objective goal attainment optimization 
% problem.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT)
%   tries to make the objective functions (F) supplied by the function FUN
%   attain the goals (GOAL) by varying X. The goals are weighted according 
%   to WEIGHT. In doing so the following nonlinear programming problem is 
%   solved:
%            min     { GAMMA :  F(X)-WEIGHT.*GAMMA<=GOAL } 
%          X,GAMMA  
%
%   FUN accepts input X and returns a vector (matrix) of function values F 
%   evaluated at X. X0 may be a scalar, vector, or matrix.  
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B) solves the goal attainment 
%   problem subject to the linear inequalities A*X <= B.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq) solves the goal
%   attainment problem subject to the linear equalities Aeq*X = Beq as
%   well.  
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB) defines a set of 
%   lower and upper bounds on the design variables, X, so that the solution
%   is in the range LB <= X <= UB. Use empty matrices for LB and U if no 
%   bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; set 
%   UB(i) = Inf if X(i) is unbounded above.
%   
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON) subjects 
%   the goal attainment problem to the constraints defined in NONLCON 
%   (usually a MATLAB file: NONLCON.m). The function NONLCON should return the 
%   vectors C and Ceq, representing the nonlinear inequalities and 
%   equalities respectively, when called with feval: 
%   [C, Ceq] = feval(NONLCON,X). FGOALATTAIN optimizes such that C(X) <= 0 
%   and Ceq(X) = 0.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
%   minimizes the with default optimization parameters replaced by values
%   in the structure OPTIONS, an argument created with the OPTIMSET
%   function.  See OPTIMSET for details. Used options are Display, TolX,
%   TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, GradConstr,
%   MaxFunEvals, MaxIter, MeritFunction, GoalsExactAchieve, Diagnostics,
%   DiffMinChange, DiffMaxChange, PlotFcns, OutputFcn, and TypicalX. Use
%   the GradObj option to specify that FUN may be called with two output
%   arguments where the second, G, is the partial derivatives of the
%   function df/dX, at the point X: [F,G] = feval(FUN,X). Use the
%   GradConstr option to specify that NONLCON may be called with four
%   output arguments: [C,Ceq,GC,GCeq] = feval(NONLCON,X) where GC is the
%   partial derivatives of the constraint vector of inequalities C an GCeq
%   is the partial derivatives of the constraint vector of equalities Ceq.
%   Use OPTIONS = [] as a place holder if no options are set.
%
%   X = FGOALATTAIN(PROBLEM) solves the goal attainment problem defined in 
%   PROBLEM. PROBLEM is a structure with the function FUN in 
%   PROBLEM.objective, the start point in PROBLEM.x0, the 'goal' vector in 
%   PROBLEM.goal, the 'weight' vector in PROBLEM.weight, the linear 
%   inequality constraints in PROBLEM.Aineq and PROBLEM.bineq, the linear 
%   equality constraints in PROBLEM.Aeq and PROBLEM.beq, the lower bounds 
%   in PROBLEM.lb, the upper bounds in PROBLEM.ub, the nonlinear constraint
%   function in PROBLEM.nonlcon, the options structure in PROBLEM.options, 
%   and solver name 'fgoalattain' in PROBLEM.solver. Use this syntax to 
%   solve at the command line a problem exported from OPTIMTOOL. The 
%   structure PROBLEM must have all the fields.
%
%   [X,FVAL] = FGOALATTAIN(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,ATTAINFACTOR] = FGOALATTAIN(FUN,X0,...) returns the attainment
%   factor at the solution X. If ATTAINFACTOR is negative, the goals have
%   been over- achieved; if ATTAINFACTOR is positive, the goals have been
%   under-achieved.
%
%   [X,FVAL,ATTAINFACTOR,EXITFLAG] = FGOALATTAIN(FUN,X0,...) returns an 
%   EXITFLAG that describes the exit condition of FGOALATTAIN. Possible 
%   values of EXITFLAG and the corresponding exit conditions are listed
%   below. See the documentation for a complete description.
%
%     1  FGOALATTAIN converged to a solution.
%     4  Computed search direction too small.
%     5  Predicted change in ATTAINFACTOR too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  No feasible point found.
%   
%   [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT] = FGOALATTAIN(FUN,X0,...) returns 
%   a structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the norm of the final step in OUTPUT.stepsize, the 
%   final line search steplength in OUTPUT.lssteplength, the algorithm used
%   in OUTPUT.algorithm, the first-order optimality in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
% 
%   [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = FGOALATTAIN(FUN,X0,...)
%   returns the Lagrange multiplier at the solution X: LAMBDA.lower for
%   LB, LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear
%   inequalities, LAMBDA.eqlin is for the linear equalities,
%   LAMBDA.ineqnonlin is for the nonlinear inequalities, and
%   LAMBDA.eqnonlin is for the nonlinear equalities.
%
%   See also OPTIMSET, OPTIMGET.

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.21 $  $Date: 2011/05/09 01:05:46 $

% ---------------------More Details---------------------------
% [x]=fgoalattain(F,x,GOAL,WEIGHT,[],[],[],[],[],[],[],OPTIONS)
% Solves the goal attainment problem where:
%
%  X  Is a set of design parameters which can be varied.
%  F  Is a set of objectives which are dependent on X.
%  GOAL Set of design goals. The optimizer will try to make 
%         F<GOAL, F=GOAL, or F>GOAL depending on the formulation.
%  WEIGHT Set of weighting parameters which determine the 
%         relative under or over achievement of the objectives.
%         Notes:
%           1.Setting WEIGHT=abs(GOAL)  will try to make the objectives
%             less than the goals resulting in roughly the same 
%             percentage under or over achievement of the goals.
%             Note: use WEIGHT 1 for GOALS that are 0 (see Note 3 below).
%           2. Setting WEIGHT=-abs(GOAL) will try to make the objectives
%              greater then the goals resulting in roughly the same percentage 
%              under- or over-achievement in the goals.
%             Note: use WEIGHT 1 for GOALS that are 0 (see Note 3 below).
%           3. Setting WEIGHT(i)=0  indicates a hard constraint.
%              i.e. F<=GOAL.
%  OPTIONS.GoalsExactAchieve indicates the number of objectives for which it is
%      required for the objectives (F) to equal the goals (GOAL). 
%      Such objectives should be partitioned into the first few 
%      elements of F.
%      The remaining parameters determine tolerance settings.
%          
%
%
defaultopt = struct( ...
    'DerivativeCheck','off', ...
    'Diagnostics','off', ...
    'DiffMaxChange',Inf, ...
    'DiffMinChange',0, ...
    'Display','final', ...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward', ...
    'FunValCheck','off', ...
    'GoalsExactAchieve',0, ...
    'GradConstr','off', ...
    'GradObj','off', ...
    'Hessian','off', ...
    'LargeScale','off', ...
    'MaxFunEvals','100*numberOfVariables', ...
    'MaxIter',400, ...
    'MaxSQPIter','10*max(numberOfVariables,numberOfInequalities+numberOfBounds)', ...
    'MeritFunction','multiobj', ...
    'NoStopIfFlatInfeas','off', ...
    'OutputFcn',[], ...
    'PhaseOneTotalScaling','off', ...
    'PlotFcns',[], ...
    'RelLineSrchBnd',[], ...
    'RelLineSrchBndDuration',1, ...
    'TolCon',1e-6, ...
    'TolConSQP',1e-6, ...
    'TolFun',1e-6, ...
    'TolX',1e-6, ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel','never' ...
    );

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 12
    options = [];
    if nargin < 11
        NONLCON = [];
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
end

% Detect problem structure input
problemInput = false;
if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fgoalattain:InputArg'));
    end
end

if nargin < 4 && ~problemInput
    error(message('optim:fgoalattain:NotEnoughInputs'))
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB);
catch ME
    if strcmp(ME.identifier,'MATLAB:datatypes:superiorfloat')
        dataType = 'notDouble';
    end
end

if ~strcmp(dataType,'double')
    error(message('optim:fgoalattain:NonDoubleInput'))
end

initVals.xOrigShape = x;
[sizes.xRows,sizes.xCols] = size(x);
xnew = [x(:); 0];
numberOfVariablesplus1 = length(xnew);
sizes.nVar = numberOfVariablesplus1 - 1;
WEIGHT = WEIGHT(:);
GOAL = GOAL(:);

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

[xnew(1:sizes.nVar),l,u,msg] = checkbounds(xnew(1:sizes.nVar),LB,UB,sizes.nVar);
if ~isempty(msg)
    EXITFLAG = -2;
    [FVAL,ATTAINFACTOR,LAMBDA] = deal([]);
    OUTPUT.iterations = 0;
    OUTPUT.funcCount = 0;
    OUTPUT.stepsize = [];
    OUTPUT.lssteplength = [];
    OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';
    OUTPUT.firstorderopt = [];
    OUTPUT.constrviolation = [];
    OUTPUT.message = msg;
    x(:) = xnew(1:sizes.nVar);
    if verbosity > 0
        disp(msg)
    end
    return
end

neqgoals = optimget(options, 'GoalsExactAchieve',defaultopt,'fast');
% flags.meritFunction is 1 unless changed by user to fmincon merit function;
% formerly options(7)
% 0 uses the fmincon single-objective merit and Hess; 1 is the default
flags.meritFunction = strcmp(optimget(options,'MeritFunction',defaultopt,'fast'),'multiobj');

lenVarIn = length(varargin);
% goalcon and goalfun also take:
% neqgoals,funfcn,gradfcn,WEIGHT,GOAL,x,errCheck
goalargs = 7; 

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
% Gather options needed for finitedifferences
% Write checked DiffMaxChange, DiffMinChage, FinDiffType, FinDiffRelStep,
% GradObj and GradConstr options back into struct for later use
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast');
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
if options.DiffMinChange >= options.DiffMaxChange
    error(message('optim:fgoalattain:DiffChangesInconsistent', sprintf( '%0.5g', options.DiffMinChange ), sprintf( '%0.5g', options.DiffMaxChange )))
end
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx(:);
options = validateFinDiffRelStep(sizes.nVar,options,defaultopt);

options.GradObj = optimget(options,'GradObj',defaultopt,'fast');
options.GradConstr = optimget(options,'GradConstr',defaultopt,'fast');

flags.grad = strcmp(options.GradObj,'on');
flags.gradconst = strcmp(options.GradConstr,'on');
if strcmpi(optimget(options,'Hessian',defaultopt,'fast'),'on')
    warning(message('optim:fgoalattain:UserHessNotUsed'))
end
flags.hess = false;

constflag = ~isempty(NONLCON);

% If nonlinear constraints exist, need either both function and constraint
% gradients, or none
if constflag
    flags.gradconst = flags.grad && flags.gradconst;
else % No user nonlinear constraints
    flags.gradconst = flags.grad;
end
flags.grad = true; % Always can compute gradient of goalfun since based on x

% Update options GradObj and GradConstr to reflect the update for the
% constraint function
if ~flags.gradconst
    options.GradObj = 'off';
    options.GradConstr = 'off';
end

% There is only one algorithm used by fgoalattain
if strcmpi(optimget(options,'LargeScale',defaultopt,'fast'),'on')
    warning(message('optim:fgoalattain:NoLargeScale'))
end

% Convert to inline function as needed
% FUN is called from goalcon; goalfun is based only on x
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    % Pass flags.gradconst as the flag which tells whether or not to
    % evaluate gradients from the user function. flags.grad is meant for
    % goalfun and is always set to true for this problem.
    funfcn = optimfcnchk(FUN,'goalcon',length(varargin),funValCheck, ...
        flags.gradconst,flags.hess);
else
    error(message('optim:fgoalattain:InvalidFUN'))
end

if constflag % NONLCON is non-empty
   confcn = optimfcnchk(NONLCON,'goalcon',length(varargin),funValCheck, ...
       flags.gradconst,false,true);
else
   confcn{1} = '';
end

% Pass in false for funValCheck argument as goalfun/goalcon is not a user function
ffun = optimfcnchk(@goalfun,'fgoalattain',lenVarIn+goalargs,false,flags.grad);
cfun = optimfcnchk(@goalcon,'fgoalattain',lenVarIn+goalargs,false,flags.gradconst,false,true); 

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
x(:) = xnew(1:end-1);
sizes.nFun = length(GOAL); % Assume the length of GOAL is same as length
                           % of user function; we will verify this later.

% Check if neqgoals (GoalsExactAchieve) is less or equal to the length of user function                           
if neqgoals > sizes.nFun
    warning(message('optim:fgoalattain:InconsistentNumEqGoal'))
    % The number of goals to be achieved exactly can be at most equal to the
    % length of user objective function.
    neqgoals = sizes.nFun;
end                         
                           
if length(WEIGHT) ~= length(GOAL)
     error(message('optim:fgoalattain:InvalidWeightAndGoalSizes'))
end

initVals.g = zeros(numberOfVariablesplus1,1);
initVals.H = [];
errCheck = true; % Perform error checking on initial function evaluations

extravarargin = {neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:}}; 
% Evaluate goal function
switch ffun{1}
    case 'fun'
        initVals.f = feval(ffun{3},xnew,extravarargin{:});
    case 'fungrad'
        [initVals.f,initVals.g] = feval(ffun{3},xnew,extravarargin{:});
    otherwise
        error(message('optim:fgoalattain:InvalidCalltype'))
end


% Evaluate goal constraints
switch cfun{1}
    case 'fun'
        [ctmp,ceqtmp] = feval(cfun{3},xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
        initVals.gnc = zeros(numberOfVariablesplus1,length(initVals.ncineq));
        initVals.gnceq = zeros(numberOfVariablesplus1,length(initVals.nceq));
    case 'fungrad'
        [ctmp,ceqtmp,initVals.gnc,initVals.gnceq] = feval(cfun{3},xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
    otherwise
        error(message('optim:fgoalattain:InvalidCalltype'))
end

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

sizes.mNonlinEq = length(initVals.nceq);
sizes.mNonlinIneq = length(initVals.ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);

if Aeqcol ~= sizes.nVar
   error(message('optim:fgoalattain:InvalidSizeOfAeq', sizes.nVar))
end
if Acol ~= sizes.nVar
   error(message('optim:fgoalattain:InvalidSizeOfA', sizes.nVar))
end

just_user_constraints = sizes.mNonlinIneq - sizes.nFun - neqgoals;
OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';  % override nlconst output

if diagnostics
    % Do diagnostics on information so far
    diagnose('fgoalattain',OUTPUT,flags.gradconst,flags.hess,constflag,flags.gradconst,...
        xnew(1:end-1),sizes.mNonlinEq,just_user_constraints,lin_eq,lin_ineq,LB,UB,funfcn,confcn);
end

% Add extra column to account for extra xnew component
A = [A,zeros(lin_ineq,1)];
Aeq = [Aeq,zeros(lin_eq,1)];

% Only need to perform error checking on initial function evaluations
errCheck = false;

% Convert function handles to anonymous functions with additional arguments
% in its workspace. Even though ffun and cfun are internal functions, put fevals
% here for consistency.
ffun{3} = @(y,varargin) feval(ffun{3},y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});
cfun{3} = @(y,varargin) feval(cfun{3},y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});

% Problem related data is passed to nlconst in problemInfo structure
problemInfo.nHardConstraints = neqgoals;
problemInfo.weight = WEIGHT;
problemInfo.goal = GOAL;

% Create default structure of flags for finitedifferences:
% This structure will (temporarily) ignore some of the features that are
% algorithm-specific (e.g. scaling and fault-tolerance) and can be turned
% on later for the main algorithm.
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward');
finDiffFlags.scaleObjConstr = false; % No scaling for now
finDiffFlags.chkFunEval = false;     % No fault-tolerance yet
finDiffFlags.chkComplexObj = false;  % No need to check for complex values
finDiffFlags.isGrad = false;         % Multi-objective
finDiffFlags.hasLBs = false(sizes.nVar,1);
finDiffFlags.hasUBs = false(sizes.nVar,1);
if ~isempty(l)
    finDiffFlags.hasLBs = isfinite(l);   % Finite lower bounds
end
if ~isempty(u)
    finDiffFlags.hasUBs = isfinite(u);   % Finite upper bounds
end

% Check first derivatives, if requested
if DerivativeCheck && ...               % User wants to check derivatives...
   (strcmpi(options.GradObj,'on') || ...                   % of either objective or ...
   strcmpi(options.GradConstr,'on') && just_user_constraints > 0) % nonlinear constraint function.
    validateFirstDerivatives(funfcn,confcn,xnew(1:sizes.nVar), ...
        l,u,options,finDiffFlags,sizes,varargin{:});
end

% Adjust nVar-length vectors used by finite-differencing for auxiliary variable
options.TypicalX = [typicalx(:); 1]; % add element for auxiliary variable
if finDiffFlags.fwdFinDiff
    options.FinDiffRelStep = [options.FinDiffRelStep; sqrt(eps)];
else
    options.FinDiffRelStep = [options.FinDiffRelStep; eps^(1/3)];
end
l = [l;-Inf];
u = [u; Inf];
finDiffFlags.hasLBs = [finDiffFlags.hasLBs; false];
finDiffFlags.hasUBs = [finDiffFlags.hasUBs; false];
finDiffFlags.isGrad = true;         % New formulation has single objective

[xnew,ATTAINFACTOR,LAMBDA,EXITFLAG,OUTPUT]=...
   nlconst(ffun,xnew,l,u,full(A),B,full(Aeq),Beq,cfun,options,defaultopt, ...
   finDiffFlags,verbosity,flags,initVals,problemInfo,varargin{:});

if ~isempty(LAMBDA)
    just_user_constraints = length(LAMBDA.ineqnonlin) - sizes.nFun - neqgoals;
    LAMBDA.ineqnonlin = LAMBDA.ineqnonlin(1:just_user_constraints);
    LAMBDA.lower = LAMBDA.lower(1:sizes.nVar);
    LAMBDA.upper = LAMBDA.upper(1:sizes.nVar);
end

OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';  % override nlconst output

% Evaluate user objective functions
x(:) = xnew(1:end-1);
FVAL = feval(funfcn{3},x,varargin{:});

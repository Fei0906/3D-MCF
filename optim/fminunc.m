function [x,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(FUN,x,options,varargin)
%FMINUNC finds a local minimum of a function of several variables.
%   X = FMINUNC(FUN,X0) starts at X0 and attempts to find a local minimizer
%   X of the function FUN. FUN accepts input X and returns a scalar
%   function value F evaluated at X. X0 can be a scalar, vector or matrix. 
%
%   X = FMINUNC(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, GradObj, HessPattern, Hessian, HessMult, HessUpdate,
%   InitialHessType, InitialHessMatrix, MaxFunEvals, MaxIter, DiffMinChange
%   and DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG,
%   PlotFcns, OutputFcn, and TypicalX. Use the GradObj option to specify
%   that FUN also returns a second output argument G that is the partial
%   derivatives of the function df/dX, at the point X. Use the Hessian
%   option to specify that FUN also returns a third output argument H that
%   is the 2nd partial derivatives of the function (the Hessian) at the
%   point X. The Hessian is only used by the large-scale algorithm. 
%
%   X = FMINUNC(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminunc' in PROBLEM.solver. Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. The structure PROBLEM 
%   must have all the fields.
%
%   [X,FVAL] = FMINUNC(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINUNC(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of FMINUNC. Possible values of EXITFLAG 
%   and the corresponding exit conditions are listed below. See the
%   documentation for a complete description.
%
%     1  Magnitude of gradient small enough. 
%     2  Change in X too small.
%     3  Change in objective function too small.
%     5  Cannot decrease function along search direction.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -3  Problem seems unbounded. 
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINUNC(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, the number of CG iterations (if used) in
%   OUTPUT.cgiterations, the first-order optimality (if used) in
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINUNC(FUN,X0,...) returns the value 
%   of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = FMINUNC(FUN,X0,...) returns the 
%   value of the Hessian of the objective function FUN at the solution X.
%
%   Examples
%     FUN can be specified using @:
%        X = fminunc(@myfun,2)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x) + 3;
%
%     To minimize this function with the gradient provided, modify
%     the function myfun so the gradient is the second output argument:
%        function [f,g] = myfun(x)
%         f = sin(x) + 3;
%         g = cos(x);
%     and indicate the gradient value is available by creating an options
%     structure with OPTIONS.GradObj set to 'on' (using OPTIMSET):
%        options = optimset('GradObj','on');
%        x = fminunc(@myfun,4,options);
%
%     FUN can also be an anonymous function:
%        x = fminunc(@(x) 5*x(1)^2 + x(2)^2,[5;1])
%
%   If FUN is parameterized, you can use anonymous functions to capture the
%   problem-dependent parameters. Suppose you want to minimize the 
%   objective given in the function myfun, which is parameterized by its 
%   second argument c. Here myfun is a MATLAB file function such as
%
%     function [f,g] = myfun(x,c)
%
%     f = c*x(1)^2 + 2*x(1)*x(2) + x(2)^2; % function
%     g = [2*c*x(1) + 2*x(2)               % gradient
%          2*x(1) + 2*x(2)];
%
%   To optimize for a specific value of c, first assign the value to c. 
%   Then create a one-argument anonymous function that captures that value 
%   of c and calls myfun with two arguments. Finally, pass this anonymous 
%   function to FMINUNC:
%
%     c = 3;                              % define parameter first
%     options = optimset('GradObj','on'); % indicate gradient is provided 
%     x = fminunc(@(x) myfun(x,c),[1;1],options)
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, @, INLINE.

%   When options.LargeScale=='on', the algorithm is a trust-region method.
%   When options.LargeScale=='off', the algorithm is the BFGS Quasi-Newton 
%   method with a mixed quadratic and cubic line search procedure. 

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.19 $  $Date: 2011/10/15 01:57:40 $
 
% ------------Initialization----------------
defaultopt = struct( ...
    'DerivativeCheck','off', ...   
    'Diagnostics','off', ...
    'DiffMaxChange',Inf, ...
    'DiffMinChange',0, ...
    'Display','final', ...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward', ...
    'FunValCheck','off', ...
    'GradObj','off', ...
    'Hessian','off', ...
    'HessMult',[], ...
    'HessPattern','sparse(ones(numberOfVariables))', ...
    'HessUpdate','bfgs', ...
    'InitialHessType','scaled-identity', ...
    'InitialHessMatrix',[], ...
    'LargeScale','on', ...
    'MaxFunEvals','100*numberOfVariables', ...
    'MaxIter',400, ...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
    'ObjectiveLimit', -1e20, ...
    'OutputFcn',[], ...
    'PlotFcns',[], ...
    'PrecondBandWidth',0, ...
    'TolFun',1e-6, ...
    'TolPCG',0.1, ...
    'TolX',1e-6, ...
    'TypicalX','ones(numberOfVariables,1)' ...
    ); 

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, options=[]; end 

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fminunc:InputArg'));
    end
end

if nargin == 0 
  error(message('optim:fminunc:NotEnoughInputs'))
end

if nargout > 5
  flags.computeHessian = true;
else
  flags.computeHessian = false;    
end

% Check for non-double inputs
if ~isa(x,'double')
  error(message('optim:fminunc:NonDoubleInput'))
end

XOUT=x(:);
sizes.nVar = length(XOUT);
sizes.mNonlinIneq = 0;
sizes.mNonlinEq = 0;
[sizes.xRows,sizes.xCols] = size(x);

medium = 'medium-scale: Quasi-Newton line search'; 
large = 'large-scale: trust-region Newton'; 

display = optimget(options,'Display',defaultopt,'fast');
flags.detailedExitMsg = ~isempty(strfind(display,'detailed'));
switch display
case {'off','none'}
   flags.verbosity = 0;
case {'notify','notify-detailed'}
   flags.verbosity = 1;  
case {'final','final-detailed'}
   flags.verbosity = 2;   
case {'iter','iter-detailed'}
   flags.verbosity = 3;
case 'testing'
   flags.verbosity = Inf;
otherwise
   flags.verbosity = 2;
end
diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');

% Check options needed for Derivative Check
options.GradObj = optimget(options,'GradObj',defaultopt,'fast');
options.GradConstr = 'off';
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');

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

DerivativeCheck = strcmpi(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
gradflag =  strcmp(options.GradObj,'on');
Hessian = optimget(options,'Hessian',defaultopt,'fast');

% line_search: 0 means trust-region, 1 means line-search
line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off');

if ( strcmpi(Hessian,'on') || strcmpi(Hessian,'user-supplied') )
    hessflag = true;
elseif strcmpi(Hessian,'off') || strcmpi(Hessian,'fin-diff-grads')
    hessflag = false;
else
    % If calling large-scale algorithm with an unavailable Hessian option value,
    % issue informative error message
    if ~line_search
        error(message('optim:fminunc:BadTRReflectHessianValue'))
    end
end

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
flags.computeLambda = 0;

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   funfcn = optimfcnchk(FUN,'fminunc',length(varargin),funValCheck,gradflag,hessflag);
else
   error(message('optim:fminunc:InvalidFUN'))
end

GRAD = zeros(sizes.nVar,1);
HESS = [];

switch funfcn{1}
case 'fun'
    try
        f = feval(funfcn{3},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:ObjectiveError', ...
            getString(message('optim:fminunc:ObjectiveError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
case 'fungrad'
    try
        [f,GRAD(:)] = feval(funfcn{3},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:ObjectiveError', ...
            getString(message('optim:fminunc:ObjectiveError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
case 'fungradhess'
    try
      [f,GRAD(:),HESS] = feval(funfcn{3},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:ObjectiveError', ...
            getString(message('optim:fminunc:ObjectiveError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
case 'fun_then_grad'
    try
        f = feval(funfcn{3},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:ObjectiveError', ...
            getString(message('optim:fminunc:ObjectiveError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
    try
        GRAD(:) = feval(funfcn{4},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:GradientError', ...
            getString(message('optim:fminunc:GradientError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
case 'fun_then_grad_then_hess'
    try
     f = feval(funfcn{3},x,varargin{:}); 
      catch userFcn_ME
        optim_ME = MException('optim:fminunc:ObjectiveError', ...
            getString(message('optim:fminunc:ObjectiveError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
     try
   GRAD(:) = feval(funfcn{4},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:GradientError', ...
            getString(message('optim:fminunc:GradientError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
   
   try
   HESS = feval(funfcn{5},x,varargin{:});
    catch userFcn_ME
        optim_ME = MException('optim:fminunc:HessianError', ...
            getString(message('optim:fminunc:HessianError')));
        userFcn_ME = addCause(userFcn_ME,optim_ME);
        rethrow(userFcn_ME)
    end
otherwise
   error(message('optim:fminunc:UndefCalltype'));
end

% Check that the objective value is a scalar
if numel(f) ~= 1
   error(message('optim:fminunc:NonScalarObj'))
end

% Determine algorithm
% If line-search and no hessian,  then call line-search algorithm
if line_search  && ...
      (~strcmpi(funfcn{1}, 'fun_then_grad_then_hess') && ~strcmpi(funfcn{1}, 'fungradhess'))
  output.algorithm = medium; 
    
  % Line-search and Hessian -- no can do, so do line-search after warning: ignoring hessian.   
elseif line_search && ...
        (strcmpi(funfcn{1}, 'fun_then_grad_then_hess') || strcmpi(funfcn{1}, 'fungradhess'))
    warning(message('optim:fminunc:HessIgnored'))
    if strcmpi(funfcn{1}, 'fun_then_grad_then_hess')
        funfcn{1} = 'fun_then_grad';
    elseif strcmpi(funfcn{1}, 'fungradhess')
        funfcn{1} = 'fungrad';
    end
    output.algorithm = medium;
    % If not line-search (trust-region) and Hessian, call trust-region   
elseif ~line_search && ...
        (strcmpi(funfcn{1}, 'fun_then_grad_then_hess') || strcmpi(funfcn{1}, 'fungradhess'))
   l=[]; u=[]; Hstr=[];
   output.algorithm = large; 
% If not line search (trust-region) and no Hessian but grad, use sparse finite-differencing.
elseif ~line_search && ...
      (strcmpi(funfcn{1}, 'fun_then_grad') || strcmpi(funfcn{1}, 'fungrad'))
   n = length(XOUT); 
   Hstr = optimget(options,'HessPattern',defaultopt,'fast');
   if ischar(Hstr) 
      if strcmpi(Hstr,'sparse(ones(numberofvariables))')
      % Put this code separate as it might generate OUT OF MEMORY error
         Hstr = sparse(ones(n));
      else
         error(message('optim:fminunc:InvalidHessPattern'))
      end
   end
   checkoptionsize('HessPattern', size(Hstr), n);
   l=[]; u=[];
   output.algorithm = large;
   
   % Trust region but no grad, no can do; warn and use line-search    
elseif ~line_search
   warning(message('optim:fminunc:SwitchingMethod'))
   output.algorithm = medium;
else
   error(message('optim:fminunc:InvalidProblem'))   
end
% Set up confcn for diagnostics and derivative check
confcn = {''};
if diagnostics
   % Do diagnostics on information so far
   constflag = false; gradconstflag = false; 
   non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
   diagnose('fminunc',output,gradflag,hessflag,constflag,gradconstflag,...
      XOUT,non_eq,non_ineq,lin_eq,lin_ineq,[],[],funfcn,confcn);
   
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
finDiffFlags.hasLBs = false(sizes.nVar,1); % No lower bounds
finDiffFlags.hasUBs = false(sizes.nVar,1); % No lower bounds

% Check derivatives
if DerivativeCheck && gradflag           % user wants to check derivatives
    validateFirstDerivatives(funfcn,confcn,XOUT,-Inf(sizes.nVar,1), ...
        Inf(sizes.nVar,1),options,finDiffFlags,sizes,varargin{:});
end

% If line-search and no hessian,  then call line-search algorithm
if strcmpi(output.algorithm, medium)
   [x,FVAL,GRAD,HESSIAN,EXITFLAG,OUTPUT] = fminusub(funfcn,x, ...
      options,defaultopt,f,GRAD,sizes,flags,finDiffFlags,varargin{:});
elseif strcmpi(output.algorithm, large)
    % Fminunc does not support output.constrviolation 
    computeConstrViolForOutput = false;
   [x,FVAL,~,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminbx(funfcn,x,l,u, ...
      flags.verbosity,options,defaultopt,flags.computeLambda,f,GRAD,HESS,Hstr, ...
      flags.detailedExitMsg,computeConstrViolForOutput,varargin{:});
   OUTPUT.algorithm = large; % override sfminbx output: not using the reflective 
                             % part of the method   
end


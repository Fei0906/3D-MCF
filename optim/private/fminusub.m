function [x,f,grad,hessian,exitflag,output] = ...
    fminusub(funfcn,x,options,defaultopt,f,grad,sizes,flags,finDiffFlags,varargin)
% FMINUSUB finds the minimizer x of a function funfcn of several variables. 
% On input, x is the initial guess, f and grad are the values of the
% function and the gradient, respectively, both evaluated at the initial
% guess x. 
% On output, x is the computed solution, f and grad are the values of the
% function and the gradient, evaluated at the computed solution x. The
% output variable hessian is the finite-difference Hessian.
 
%   Copyright 1990-2011 The MathWorks, Inc. 
%   $Revision: 1.1.6.26 $  $Date: 2011/10/15 01:57:43 $

% Check to see if the objective value at the initial point is 
% well-defined.  If not, then terminate immediately.
if ~isfinite(f) || ~isreal(f)
    error(message('optim:fminusub:UsrObjUndefAtX0','Fminunc'))
end
% If user-supplied derivatives at x0 are not defined, terminate immediately, as there is
% nothing we can do to move away from x0. 
if any(~isfinite(grad)) || any(~isreal(grad))
    error(message('optim:fminusub:GradUndefAtX0','Fminunc'))
end

%
% Initialization
%
verbosity = flags.verbosity;
detailedExitMsg = flags.detailedExitMsg;

x = x(:);                % Reshape x to a column vector
numberOfVariables = length(x);
initialHessIsAScalar = [];
exitflagLnSrch = [];     % define in case x0 is solution and lineSearch never called
dir = [];                % define for last call to outputFcn in case x0 is solution
formatstr = ' %5.0f       %5.0f    %13.6g  %13.6g   %12.3g  %s\n';

% Line search parameters: rho < 1/2 and rho < sigma < 1. Typical values are
% rho = 0.01 and sigma = 0.9.
rho = 0.01; sigma = 0.9;

% Initialize fault tolerance structure
faultTolStruct.undefObj = false;             % ran into an undefined objective value?
faultTolStruct.currTrialWellDefined = true;  % current obj value well defined?
faultTolStruct.undefObjValue = [];           % actual undefined objective value (last one if there were many)
faultTolStruct.relStepSize = [];             % relative step size

% Read in options
gradflag =  strcmp(optimget(options,'GradObj',defaultopt,'fast'),'on');
TolX = optimget(options,'TolX',defaultopt,'fast');

HessUpdate = optimget(options,'HessUpdate',defaultopt,'fast'); 
InitialHessType = optimget(options,'InitialHessType',defaultopt,'fast');
InitialHessMatrix = optimget(options,'InitialHessMatrix',defaultopt,'fast');

fminimum = optimget(options,'ObjectiveLimit',defaultopt,'fast'); 

if strcmpi(InitialHessType,'user-supplied')
    if isempty(InitialHessMatrix)    
        warning(message('optim:fminusub:ResettingToInitialHessType'))
        InitialHessType = 'identity';    
    else
        % Determine size of InitialHessMatrix
        [ihRows,ihCols] = size(InitialHessMatrix);
        if (ihRows==numberOfVariables && ihCols==1) || (ihRows==1 && ihCols==numberOfVariables)
          initialHessIsAScalar = false;
        elseif ihRows==1 && ihCols==1
          initialHessIsAScalar = true;          
        else
          error(message('optim:fminusub:InitialHessMatrixSize'))
        end    
    end
end

TolFun = optimget(options,'TolFun',defaultopt,'fast');
maxFunEvals = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxIter = optimget(options,'MaxIter',defaultopt,'fast');

if ischar(maxFunEvals)
    if isequal(lower(maxFunEvals),'100*numberofvariables')
        maxFunEvals = 100*numberOfVariables;
    else
        error(message('optim:fminusub:MaxFunEvalsIntOrDefault'))
    end
end

% Update structure of flags for finitedifferences
finDiffFlags.chkFunEval = true; % Validate function values
finDiffFlags.chkComplexObj = true; % Check whether objective function values are complex if chkFunEval is true

% Prepare strings to give feedback to users on options they have or have not set.
% These are used in the exit messages.
optionFeedback = createOptionFeedback(options);

% Output function
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = reshape(x,sizes.xRows,sizes.xCols); % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for
    % OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Plot functions
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = reshape(x,sizes.xRows,sizes.xCols); % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for
    % PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

funcCount = 1; % function evaluated in FMINUNC
iter = 0;

% Initialize output alpha: if x0 is the solution, alpha = [] is returned in
% output structure
alpha = []; 
hessUpdateMsg = [];       

% Compute finite difference gradient at initial point, if needed
if ~gradflag
    grad = zeros(numberOfVariables,1); % pre-allocate finite-difference gradient
    [grad,~,~,numEvals,evalOKFinDiff] = finitedifferences(x,funfcn{3},[],[],[],f,[],[], ...
        1:numberOfVariables,options,sizes,grad,[],[],finDiffFlags,[],varargin{:});
    
    % If finite difference gradient at x0 is not defined, terminate immediately,
    % as there is nothing the algorithm can do to move away from x0.
    if ~evalOKFinDiff
        error(message('optim:fminusub:DerivUndefAtX0','Fminunc'))
    end

    funcCount = funcCount + numEvals;
end

% Norm of initial gradient, used in stopping tests
g0Norm = norm(grad,Inf); 

% Initialize the output function.
if haveoutputfcn || haveplotfcn
  [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'init',iter,funcCount, ...
        f,[],grad,[],varargin{:});
  if stop
    [x,f,exitflag,output,grad,hessian] = cleanUpInterrupt(xOutputfcn,optimValues,verbosity,detailedExitMsg);
    return;
  end
end

% Print output header
if verbosity > 2
  fprintf(['                                                        First-order \n',...
  ' Iteration  Func-count       f(x)        Step-size       optimality\n']);
end

% Display 0th iteration quantities
if verbosity > 2
  fprintf(' %5.0f       %5.0f    %13.6g                  %12.3g\n',iter,funcCount,f,g0Norm);
end

% OutputFcn call 0th iteration
if haveoutputfcn || haveplotfcn
  [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'iter',iter,funcCount, ...
        f,[],grad,[],varargin{:});
  if stop  % Stop per user request.
    [x,f,exitflag,output,grad,hessian] = cleanUpInterrupt(xOutputfcn,optimValues,verbosity,detailedExitMsg);
    return;
  end
end

% Check convergence at initial point
[done,exitflag,outMessage] = initialTestStop(f,g0Norm,TolFun,fminimum, ...
    detailedExitMsg,verbosity,optionFeedback); 

% Form initial inverse Hessian approximation
if ~done                  
  H = initialQuasiNewtonMatrix(InitialHessType,InitialHessMatrix, ...
                    HessUpdate,initialHessIsAScalar,sizes.nVar);
end                       
%                     
% Main loop
%
while ~done
    iter = iter + 1;
    
    % Form search direction
    dir = -H*grad;
    dirDerivative = grad'*dir; 
    
    % Perform line search along dir
    alpha1 = 1;
    if iter == 1 
        alpha1 = min(1/g0Norm,1); 
    end  
    fOld = f; gradOld = grad; alphaOld = alpha;

    % During line search, don't exceed the overall total maxFunEvals.
    maxFunEvalsLnSrch = maxFunEvals - funcCount;

    [alpha,f,grad,exitflagLnSrch,funcCountLnSrch,faultTolStruct] = ... 
          lineSearch(funfcn,x,dir,f,dirDerivative, ...
          alpha1,rho,sigma,fminimum,maxFunEvalsLnSrch,eps(max(1,abs(f))), ...
          options,finDiffFlags,sizes,grad,TolX,faultTolStruct,varargin{:});
    funcCount = funcCount + funcCountLnSrch;
    
    % If undefined objective value encountered, display message
    if verbosity > 2
        if faultTolStruct.undefObj
            fprintf(getString(message('optimlib:commonMsgs:ObjInfNaNComplex', ...
                formatUndefValue(faultTolStruct.undefObjValue))));
        end
    end
    
    % If line search was not able to find point that satisfies both Wolfe
    % conditions, stop main iteration. If the trial point returned by unsuccessful
    % line search call is no better than previous one, restore to previoius
    % trial point.
    fvalNoBetter = ~(faultTolStruct.currTrialWellDefined) || f >= fOld;
    if exitflagLnSrch < 0 && fvalNoBetter
      % Restore previous values
      alpha = alphaOld;
      f = fOld;
      grad = gradOld;
      break
    end
    
    % Update iterate
    deltaX = alpha*dir;
    x = x + deltaX;
    
    % Display iteration quantities
    if verbosity > 2
      % Print header periodically
      if mod(iter,20) == 0
        fprintf(['                                                        First-order \n', ...
            ' Iteration  Func-count       f(x)        Step-size       optimality\n']);        
      end
        fprintf(formatstr,iter,funcCount,f,alpha,norm(grad,inf),hessUpdateMsg)
    end

    % OutputFcn call
    if haveoutputfcn || haveplotfcn
      [xOutputfcn,optimValues,stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'iter',iter,funcCount, ...
           f,alpha,grad,dir,varargin{:});      
      if stop  % Stop per user request.
        [x,f,exitflag,output,grad,hessian] = ...
                cleanUpInterrupt(xOutputfcn,optimValues,verbosity,detailedExitMsg);
        return;
      end
    end
        
    [done,exitflag,outMessage] = testStop(x,f,deltaX,iter,funcCount,TolX,TolFun,maxIter, ...
                                     maxFunEvals,fminimum,grad,g0Norm,exitflagLnSrch, ...
                                     detailedExitMsg,verbosity,optionFeedback);
    if ~done
        % Update quasi-Newton matrix.
        [H,hessUpdateMsg] = updateQuasiNewtonMatrix(H,deltaX,grad-gradOld,HessUpdate, ...
                                            InitialHessType,iter);
    end
    
end % of while

% Handle cases in which line search didn't terminate normally
if funcCount >= maxFunEvals
  exitflag = 0;
  outMessage = createExitMsg('fminusub',exitflag,verbosity > 0,detailedExitMsg, ...
      'fminunc',[],optionFeedback.MaxFunEvals,maxFunEvals);
elseif exitflagLnSrch == -2
  exitflag = 5;
  outMessage = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
      'fminunc');
elseif exitflagLnSrch == -3
  exitflag = 2;
  outMessage = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
      'fminunc',faultTolStruct.relStepSize,optionFeedback.TolX,TolX);
end
  
% Compute finite-difference Hessian only if asked for in output
if flags.computeHessian
  if verbosity > 1
      if ~gradflag
        fprintf(getString(message('optim:fminusub:FinDiffHessMsg')));          
        % If problem large, estimating the finite difference Hessian with
        % only function values may take time
        if sizes.nVar >= 100
          fprintf(getString(message('optim:fminusub:LargeFinDiffHessMsg')));            
        end
      end
  end
  hessian = finDiffHessian(funfcn,x,sizes.xRows,sizes.xCols,sizes.nVar,gradflag,f, ...
                           grad,options,varargin{:});
else
  hessian = [];
end  

% OutputFcn call
if haveoutputfcn || haveplotfcn
  callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',iter,funcCount, ...
      f,alpha,grad,dir,varargin{:});      
end   

x = reshape(x,sizes.xRows,sizes.xCols); % restore user shape
output.iterations = iter;
output.funcCount = funcCount;
output.stepsize = alpha;
output.firstorderopt = norm(grad,inf);
output.algorithm = 'medium-scale: Quasi-Newton line search';
output.message = outMessage;
%--------------------------------------------------------------------------
function [done,exitflag,msg] = testStop(x,fval,deltaX,iter,funcCount,TolX, ...
                        TolFun,maxIter,maxFunEvals,fminimum,grad,g0Norm, ...
                        exitflagLnSrch,detailedExitMsg,verbosity,optionFeedback)
%
% testStop checks if the stopping conditions are met

if norm(grad,Inf) < TolFun*(1+g0Norm)
     done = true;
     exitflag = 1;
     msg = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
         'fminunc',norm(grad,Inf)/(1+g0Norm),optionFeedback.TolFun,TolFun);
elseif exitflagLnSrch == 1 % Problem unbounded
    done = true;
    exitflag = -3;
    msg = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
         'fminunc',fval,optionFeedback.ObjectiveLimit,fminimum);     
elseif relDeltaXFminunc(x,deltaX) < TolX 
    done = true;
    exitflag = 2;
    msg = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
        'fminunc',relDeltaXFminunc(x,deltaX),optionFeedback.TolX,TolX);
elseif exitflagLnSrch == -2 % Line search could not reduce function value any more
    done = true;    % Exit Message will be handled outside of main loop
    exitflag = 5;
    msg = '';
elseif funcCount >= maxFunEvals
     done = true;   % Exit Message will be handled outside of main loop
     exitflag = 0;
     msg = '';
elseif iter > maxIter 
     done = true;
     exitflag = 0;
     % Call createExitMsg with createExitMsgExitflag = 10 for MaxIter exceeded
     msg = createExitMsg('fminusub',10,verbosity > 0,detailedExitMsg, ...
         'fminunc',[],optionFeedback.MaxIter,maxIter);
else
   exitflag = [];
   done = false;
   msg = [];
end

%-------------------------------------------------------------------------------
function [done,exitflag,msg] = initialTestStop(fval,g0Norm,TolFun,fminimum, ...
    detailedExitMsg,verbosity,optionFeedback)
%
% initialTestStop checks if the starting point satisfies a convergence
% criterion.

if g0Norm < TolFun
  % Call createExitMsg with exitflag = 100 for an optimal x0
  msg = createExitMsg('fminusub',100,verbosity > 1,detailedExitMsg, ...
      'fminunc',g0Norm,optionFeedback.TolFun,TolFun);
  done = true;
  exitflag = 1;
elseif fval <= fminimum
  done = true;
  exitflag = -3;
  msg = createExitMsg('fminusub',exitflag,verbosity > 1,detailedExitMsg, ...
         'fminunc',fval,optionFeedback.ObjectiveLimit,fminimum);
else
  exitflag = [];
  done = false;
  msg = [];
end   

%--------------------------------------------------------------------------
function H = initialQuasiNewtonMatrix(InitialHessType,InitialHessMatrix, ...
                           HessUpdate,initialHessIsAScalar,numberOfVariables)
%initialQuasiNewtonMatrix sets the initial quasi-Newton matrix that approximates
% the inverse to the Hessian.

% Unless running steepest-descent, compute initial H
if ~strcmpi(HessUpdate,'steepdesc')                         % not steepest-descent
    if isequal(lower(InitialHessType),'identity') || ...
        isequal(lower(InitialHessType),'scaled-identity')
        % Built-in initial quasi-Newton matrix. In 'scaled-identity' case,
        % the scaling occurs right after the end of the 1st iteration.
        H = eye(numberOfVariables);
    else
        % User-supplied initial approximation to the Hessian. We invert 
        % this initial matrix because we maintain an approximation H to the
        % inverse of the Hessian. The check for InitialHessMatrix > 0 was
        % already done in optimset.m
        if initialHessIsAScalar
            % InitialHessMatrix is a scalar
            H = 1/InitialHessMatrix*eye(numberOfVariables);
        else
            % InitialHessMatrix is a vector
            H = diag(1./InitialHessMatrix);
        end
    end
else
    % Steepest-descent: H is always the identity
    H = eye(numberOfVariables);
end

%--------------------------------------------------------------------------
function [H,msg] = updateQuasiNewtonMatrix(H,deltaX,deltaGrad,HessUpdate, ...
                                           InitialHessType,iter)
% UPDATEQUASINEWTONMATRIX updates the quasi-Newton matrix that approximates
% the inverse to the Hessian.

deltaXDeltaGrad = deltaX'*deltaGrad;
updateOk = deltaXDeltaGrad >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaGrad) );
if iter == 1 && strcmpi(InitialHessType,'scaled-identity')
    if updateOk
        % Reset the initial quasi-Newton matrix to a scaled identity aimed
        % at reflecting the size of the inverse true Hessian
        H = deltaXDeltaGrad/(deltaGrad'*deltaGrad)*eye(length(deltaX));
    end
end

if strcmpi(HessUpdate,'bfgs')
  if updateOk
    HdeltaGrad = H*deltaGrad;
    % BFGS update
    H = H + (1 + deltaGrad'*HdeltaGrad/deltaXDeltaGrad) * ...
        deltaX*deltaX'/deltaXDeltaGrad - (deltaX*HdeltaGrad' + ... 
        HdeltaGrad*deltaX')/deltaXDeltaGrad;
    msg = '';
  else
    msg = 'skipped update';
  end
elseif strcmpi(HessUpdate,'dfp')
  if updateOk
    HdeltaGrad = H*deltaGrad;
    % DFP update
    H = H + deltaX*deltaX'/deltaXDeltaGrad - HdeltaGrad*HdeltaGrad'/(deltaGrad'*HdeltaGrad);    
    msg = '';
  else
    msg = 'skipped update';
  end  
elseif strcmpi(HessUpdate,'steepdesc')
  % Steepest descent
  H = eye(length(deltaX));
  msg = '';  
else
    ME = MException('optim:fminusub:UnknownHessUpdate', ...
        getString(message('optim:fminusub:UnknownHessUpdate')));                    
    throwAsCaller(ME)
end

%--------------------------------------------------------------------------  
function [Hessian,functionCalls] = finDiffHessian(funfcn,x,xRows,xCols, ...
            numberOfVariables,useGrad,f,grad,options,varargin) 
% FINDIFFHESSIAN calculates the numerical Hessian of funfcn evaluated at x
% using finite differences. 

Hessian = zeros(numberOfVariables);

if useGrad
  % Define stepsize 
  CHG = sqrt(eps)*sign(x).*max(abs(x),1); 

  % Make sure step size lies within DiffMinChange and DiffMaxChange
  CHG = sign(CHG+eps).*min(max(abs(CHG),options.DiffMinChange),options.DiffMaxChange);
  % Calculate finite difference Hessian by columns with the user-supplied
  % gradient. We use the forward difference formula.
  %
  % Hessian(:,j) = 1/h(j) * [grad(x+h(j)*ej) - grad(x)]               (1)
  for j = 1:numberOfVariables
    xplus = x; 
    xplus(j) = x(j) + CHG(j);
    % evaluate gradPlus 
    switch funfcn{1}          
     case 'fun'
      error(message('optim:fminusub:WrongUseGrad'))
     case 'fungrad'
      [~,gradPlus] = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
      gradPlus = gradPlus(:);
     case 'fun_then_grad'
      gradPlus = feval(funfcn{4},reshape(xplus,xRows,xCols),varargin{:});
      gradPlus = gradPlus(:);
     otherwise
      error(message('optim:fminusub:UndefCallType'))
    end    
    % Calculate jth column of Hessian
    Hessian(:,j) = (gradPlus - grad) / CHG(j);
  end
  % Symmetrize the Hessian
  Hessian = 0.5*(Hessian + Hessian');
  
  functionCalls = numberOfVariables;
else % of 'if useGrad'
  % Define stepsize  
  CHG = eps^(1/4)*sign(x).*max(abs(x),1);  
  
  % Make sure step size lies within DiffMinChange and DiffMaxChange
  CHG = sign(CHG+eps).*min(max(abs(CHG),options.DiffMinChange),options.DiffMaxChange);
  % Calculate the upper triangle of the finite difference Hessian element 
  % by element, using only function values. The forward difference formula 
  % we use is
  %
  % Hessian(i,j) = 1/(h(i)*h(j)) * [f(x+h(i)*ei+h(j)*ej) - f(x+h(i)*ei) 
  %                          - f(x+h(j)*ej) + f(x)]                   (2) 
  % 
  % The 3rd term in (2) is common within each column of Hessian and thus
  % can be reused. We first calculate that term for each column and store
  % it in the row vector fplus_array.
  fplus_array = zeros(1,numberOfVariables);
  for j = 1:numberOfVariables
    xplus = x;
    xplus(j) = x(j) + CHG(j);
    % evaluate  
    switch funfcn{1}
     case 'fun'
      fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});       
     case 'fungrad'
      fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
     case 'fun_then_grad'  
      fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
     otherwise
      error(message('optim:fminusub:UndefCallType'))
    end    
    fplus_array(j) = fplus;
  end
  
  for i = 1:numberOfVariables
    % For each row, calculate the 2nd term in (4). This term is common to
    % the whole row and thus it can be reused within the current row: we
    % store it in fplus_i.
    xplus = x;
    xplus(i) = x(i) + CHG(i);
    % evaluate  
    switch funfcn{1}
     case 'fun'
      fplus_i = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});        
     case 'fungrad'
      fplus_i = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
     case 'fun_then_grad'  
      fplus_i = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
     otherwise
      error(message('optim:fminusub:UndefCallType'))
    end     
 
    for j = i:numberOfVariables   % start from i: only upper triangle
      % Calculate the 1st term in (2); this term is unique for each element
      % of Hessian and thus it cannot be reused.
      xplus = x;
      xplus(i) = x(i) + CHG(i);
      xplus(j) = xplus(j) + CHG(j);
      % evaluate  
      switch funfcn{1}
       case 'fun'
        fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});        
       case 'fungrad'
        fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
       case 'fun_then_grad'  
        fplus = feval(funfcn{3},reshape(xplus,xRows,xCols),varargin{:});
       otherwise
        error(message('optim:fminusub:UndefCallType'))
      end    
      Hessian(i,j) = (fplus - fplus_i - fplus_array(j) + f)/(CHG(i)*CHG(j)); 
    end 
  end % of "for i = 1:numberOfVariables"
  % Fill in the lower triangle of the Hessian
  Hessian = Hessian + triu(Hessian,1)';
  functionCalls = 2*numberOfVariables + ...        % 2nd and 3rd terms,
      numberOfVariables*(numberOfVariables + 1)/2; % 1st term in (2)
end % of 'if useGrad'

%--------------------------------------------------------------------------
function undefFvalDisp = formatUndefValue(undefFval)
% 
    if isnan(undefFval)
        undefFvalDisp = 'NaN';
    elseif ~isfinite(undefFval)
        undefFvalDisp = 'Inf';
    else % Result is complex
        undefFvalDisp = 'complex';
    end

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,funcCount, ...
    f,alpha,grad,dir,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then
% calls the outputfcn/plotfcns.  

% state - can have the values 'init', 'iter', or 'done'. 
% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = funcCount;
optimValues.fval = f;
optimValues.stepsize = alpha;
if ~isempty(dir)
    optimValues.directionalderivative = dir'*grad;
else
    optimValues.directionalderivative = [];
end
optimValues.gradient = grad;
optimValues.searchdirection = dir;
optimValues.firstorderopt = norm(grad,Inf);
optimValues.procedure = '';
xOutputfcn(:) = x;  % Set x to have user expected size

stop = false;
% Call output function
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        % case 'interrupt' No 'interrupt' case in fminusub            
        otherwise
            error(message('optim:fminusub:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('optim:fminusub:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end
%--------------------------------------------------------------------------
function [x,fval,exitflag,output,gradient,hessian] = cleanUpInterrupt(xOutputfcn,optimValues,verbosity,detailedExitMsg)
% CLEANUPINTERRUPT updates or sets all the output arguments of NLCONST when the optimization 
% is interrupted.  The HESSIAN and LAMBDA are set to [] as they may be in a state that is 
% inconsistent with the other values since we are interrupting mid-iteration.

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn;
fval = optimValues.fval;
exitflag = -1; 
output.iterations = optimValues.iteration;
output.funcCount = optimValues.funccount;
output.stepsize = optimValues.stepsize;
output.algorithm = 'medium-scale: Quasi-Newton line search';
output.firstorderopt = optimValues.firstorderopt; 
output.cgiterations = [];
output.message = createExitMsg('fminusub',exitflag,verbosity > 0,detailedExitMsg,'fminunc');
gradient = optimValues.gradient;
hessian = []; % May be in an inconsistent state






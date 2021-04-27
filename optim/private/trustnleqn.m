function [x,Fvec,JAC,EXITFLAG,OUTPUT,msgData]= trustnleqn(funfcn,x,verbosity,gradflag, ...
  options,defaultopt,Fvec,JAC,detailedExitMsg,optionFeedback,finDiffFlags,sizes,varargin)
%TRUSTNLEQN Trust-region dogleg nonlinear systems of equation solver.
%
%   TRUSTNLEQN solves a system of nonlinear equations using a dogleg trust
%   region approach.  The algorithm implemented is similar in nature
%   to the FORTRAN program HYBRD1 of J.J. More', B.S. Garbow and K.E. 
%   Hillstrom, User Guide for MINPACK 1, Argonne National Laboratory, 
%   Rept. ANL-80-74, 1980, which itself was based on the program CALFUN 
%   of M.J.D. Powell, A Fortran subroutine for solving systems of
%   nonlinear algebraic equations, Chap. 7 in P. Rabinowitz, ed.,
%   Numerical Methods for Nonlinear Algebraic Equations, Gordon and
%   Breach, New York, 1970.

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.20 $  $Date: 2011/08/29 20:34:18 $
%
% NOTE: 'x' passed in and returned in matrix form.
%       'Fvec' passed in and returned in vector form.
%
% Throughout this routine 'x' and 'F' are matrices while
% 'xvec', 'xTrial', 'Fvec' and 'FTrial' are vectors. 
% This was done for compatibility with the 'fsolve.m' interface.

% Check to see if the function vector at the initial point is well-defined.
% If not, then terminate immediately.
if any(~isfinite(Fvec)) 
    error(message('optim:trustnleqn:UsrObjUndefAtX0'));    
end

% Define some sizes.
xvec = x(:);         % vector representation of x
% Convert values to full to avoid unnecessary sparse operation overhead
Fvec = full(Fvec); 

% Get user-defined options.
[maxfunc,maxit,tolf,tolx,giventypx,outputfcn,plotfcns] = ...
    getOpts(sizes.nVar,options,defaultopt);

typx = options.TypicalX;
scale = giventypx;  % scaling featured only enabled when typx values provided


% Handle the output function.
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot function.
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

% Initialize local arrays.
d       = zeros(sizes.nVar,1);
scalMat = ones(sizes.nVar,1); 

% Initialize some trust region parameters.
Delta    = 1e0;
DeltaMax = 1e10;
eta1     = 0.05;
eta2     = 0.9;
alpha1   = 2.5;
alpha2   = 0.25;

% Other initializations.
iter = 0;
numFevals = 1;   % computed in fsolve.m
if gradflag
  numJevals = 1; % computed in fsolve.m
else
  numJevals = 0;
end 
stepAccept = true;
normd = 0.0e0;
scalemin = eps;
scalemax = 1/scalemin;
objold = 1.0e0;
obj = 0.5*Fvec'*Fvec;  % Initial Fvec computed in fsolve.m

numFDfevals = 0;
% Compute initial finite difference Jacobian, objective and gradient.
if ~gradflag
    % User is not specifying the Jacobian. Compute finite differences of
    % full Jacobian at initial point.
    JACfindiff = zeros(sizes.nFun,sizes.nVar); % pre-allocate derivative array
    [JACfindiff,~,~,numFDfevals] = finitedifferences(x,funfcn{3},[],[],[],Fvec,[],[], ...
        1:sizes.nVar,options,sizes,JACfindiff,[],[],finDiffFlags,[],varargin{:});    
    % Set initial Jacobian to be the finite difference approximation.
    JAC = JACfindiff;
end

% Increment the number of function evaluations with those used in the
% finite difference calls.
numFevals = numFevals + numFDfevals;

grad = JAC'*Fvec;
normgradinf = norm(grad,inf);

% Print header.
header = sprintf(['\n                                         Norm of      First-order   Trust-region\n',...
                    ' Iteration  Func-count     f(x)          step         optimality    radius']);
formatstr = ' %5.0f      %5.0f   %13.6g  %13.6g   %12.3g    %12.3g';
if verbosity > 1
  disp(header);
end

% Initialize the output function.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'init',iter, ...
        numFevals,Fvec,[],[],[],Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
        msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
        return;
    end
end

% Compute initial diagonal scaling matrix.
if scale
  if giventypx && ~isempty(typx) % scale based on typx values
    typx(typx==0) = 1; % replace any zero entries with ones
    scalMat = 1./abs(typx);
  else         % scale based on norm of the Jacobian (not currently active)  
    scalMat = getscalMat(sizes.nVar,JAC,scalemin,scalemax);
  end
end

% Display initial iteration information.
formatstr0 = ' %5.0f      %5.0f   %13.6g                  %12.3g    %12.3g';
% obj is 0.5*F'*F but want to display F'*F
iterOutput0 = sprintf(formatstr0,iter,numFevals,2*obj,normgradinf,Delta);
if verbosity > 1
   disp(iterOutput0);
end
% OutputFcn call.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
        msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
        return;
    end
end


% Test convergence at initial point.
[done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolx,...
     stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
     obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity);

% Create the fault tolerance structure
faultTolStruct = createFaultTolStruct(false);
newFaultTolStruct = faultTolStruct;

% Beginning of main iteration loop.
while ~done
  iter = iter + 1;

  % Compute step, d, using dogleg approach.
  [d,quadObj,normd,normdscal] = ...
       dogleg(sizes.nVar,Fvec,JAC,grad,Delta,scalMat,varargin);

  % Compute the model reduction given by d (pred).
  pred = -quadObj;

  % Compute the trial point, xTrial.
  xTrial = xvec + d;

  % Evaluate nonlinear equations and objective at trial point.
  switch funfcn{1}
  case 'fun'
    F = feval(funfcn{3},reshape(xTrial,sizes.xRows,sizes.xCols),varargin{:});
  case 'fungrad'
    [F,JACTrial] = feval(funfcn{3},reshape(xTrial,sizes.xRows,sizes.xCols),varargin{:});
    numJevals = numJevals + 1;
  case 'fun_then_grad'
    F = feval(funfcn{3},reshape(xTrial,sizes.xRows,sizes.xCols),varargin{:}); 
  otherwise
    error(message('optim:trustnleqn:UndefinedCalltype'))
  end  
  numFevals = numFevals + 1;
  FTrial = full(F(:)); % make FTrial a vector, convert to full
  objTrial = 0.5*FTrial'*FTrial; 

  % Compute the actual reduction given by xTrial (ared).
  ared = obj - objTrial;

  % Compute ratio = ared/pred.
  if pred <= 0 % reject step
    ratio = 0;
  else
    ratio = ared/pred;
  end
  
  % Update fault tolerance structure.
  faultTolStruct = updateFaultTolStruct(faultTolStruct, objTrial, ...
      verbosity > 1);
  
  if haveoutputfcn % Call output functions (we don't call plot functions with 'interrupt' flag)
      [~, ~, stop] = callOutputAndPlotFcns(outputfcn,{},xvec,xOutputfcn,'interrupt',iter, ...
          numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
      if stop  % Stop per user request.
          [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
          msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
          return;
      end
  end
  
  if ratio > eta1 && faultTolStruct.currTrialWellDefined % accept step.

    xvec = xTrial; Fvec = FTrial; objold = obj; obj = objTrial;
    x(:) = xvec; % update matrix representation
    % Compute JAC at new point. (already computed with F if 'fungrad')

    % Compute finite difference Jacobian if needed.
    if ~gradflag
        [JACfindiff,~,~,numFDfevals] = finitedifferences(x,funfcn{3},...
            [],[],[],Fvec,[],[],1:sizes.nVar,options,sizes,JACfindiff,...
            [],[],finDiffFlags,[],varargin{:});
        numFevals = numFevals + numFDfevals;
    end

    switch funfcn{1}
        case 'fun'
            JAC = JACfindiff;
        case 'fungrad'
            JAC = JACTrial;
        case 'fun_then_grad'
            JAC = feval(funfcn{4},x,varargin{:});
            numJevals = numJevals + 1;
        otherwise
            error(message('optim:trustnleqn:UndefinedCalltype'))
    end
      
    grad = JAC'*Fvec;
    normgradinf = norm(grad,inf);

    % Update internal diagonal scaling matrix (dynamic scaling).
    if scale && ~giventypx
      scalMat = getscalMat(sizes.nVar,JAC,scalemin,scalemax);
    end

    stepAccept = true;
  else % reject step.
    stepAccept = false;
  end 

  % Print iteration statistics.
  if verbosity > 1
      if faultTolStruct.undefObj
          fprintf(getString(message('optimlib:commonMsgs:ObjInfNaNComplex', ...
              faultTolStruct.undefValue)));
      end
      % obj is 0.5*F'*F but want to display F'*F
      iterOutput = sprintf(formatstr,iter,numFevals,2*obj,normd,normgradinf,Delta);
      disp(iterOutput);
  end
  % OutputFcn call.
  if haveoutputfcn || haveplotfcn
      [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
      if stop
          [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
          msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
          return;
      end
  end

  % Update trust region radius.
  Delta = updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                      alpha1,alpha2,DeltaMax,faultTolStruct);

  % Check for termination.
  [done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolx,...
       stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
       obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity);
   
  % As the iteration has been completed, the fault tolerance
  % structure needs to be reset.
  faultTolStruct = newFaultTolStruct;
   
end

if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'done',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    % Optimization done, so ignore "stop"
end


% Optimization is finished.

% Assign output statistics.
OUTPUT.iterations = iter;
OUTPUT.funcCount = numFevals;
OUTPUT.algorithm = 'trust-region-dogleg';
OUTPUT.firstorderopt = normgradinf;

% TRUSTNLEQN finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxfunc,maxit,tolf,tolx,giventypx,outputfcn,plotfcns] = ...
                getOpts(nvar,options,defaultopt)
%getOpts gets the user-defined options for TRUSTNLEQN.

% Both Medium and Large-Scale options.
maxfunc = optimget(options,'MaxFunEvals',defaultopt,'fast');
if ischar(maxfunc)
  if isequal(lower(maxfunc),'100*numberofvariables')
    maxfunc = 100*nvar;
  else
    error(message('optim:trustnleqn:InvalidMaxFunEvals'))
  end
end
maxit = optimget(options,'MaxIter',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');

% Check if TypicalX is the default (ones) 
giventypx = any(options.TypicalX ~= 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolx,...
     stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
     obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity)
%testStop checks the termination criteria for TRUSTNLEQN.

done = false;
EXITFLAG = 0;
msgData = {};

% Check termination criteria.
if stepAccept && normgradinf < tolf
  done = true;
  EXITFLAG = 1;
  if iter == 0
      msgFlag = 100;
  else
      msgFlag = EXITFLAG;
  end
  % Setup input parameters for createExitMsg with msgFlag = 100 if x0 is
  % optimal, otherwise msgFlag = 1
  msgData = {'trustnleqn',msgFlag,verbosity > 0,detailedExitMsg,'fsolve', ...
      normgradinf,optionFeedback.TolFun,tolf,2*obj,optionFeedback.TolFun,sqrt(tolf)};
elseif iter > 1 && max(abs(d)./(abs(xvec)+1)) < max(tolx^2,eps)
   % Assign msgFlag, a unique internal exitflag, to 2 or -22 for this
   % stopping test depending on whether the result appears to be a root or
   % not.
   if 2*obj < sqrt(tolf) % fval'*fval < sqrt(tolf)
      EXITFLAG = 2; msgFlag = 2;
      dispMsg = verbosity > 0;
   else
      EXITFLAG = -2; msgFlag = -22;
      dispMsg = verbosity > 0;
   end
   % Setup input parameters for createExitMsg
   msgData = {'trustnleqn',msgFlag,dispMsg,detailedExitMsg,'fsolve', ...
       max(abs(d)./(abs(xvec)+1)),optionFeedback.TolX,max(tolx^2,eps), ...
       2*obj,optionFeedback.TolFun,sqrt(tolf)};
   done = true;
elseif iter > 1 && stepAccept && normd < 0.9*Delta ...
                && abs(objold-obj) < max(tolf^2,eps)*(1+abs(objold))
  % Assign msgFlag, a unique internal exitflag, to 3 or -23 for this
  % stopping test depending on whether the result appears to be a root or
  % not.
  if 2*obj < sqrt(tolf) % fval'*fval < sqrt(tolf)
     EXITFLAG = 3; msgFlag = 3;
     dispMsg = verbosity > 0;
  else
     EXITFLAG = -2; msgFlag = -23;
     dispMsg = verbosity > 0;
  end
  % Setup input parameters for createExitMsg
  msgData = {'trustnleqn',msgFlag,dispMsg,detailedExitMsg,'fsolve', ...
      abs(objold-obj)./(abs(objold)+1),optionFeedback.TolFun,max(tolf^2,eps), ...
      2*obj,optionFeedback.TolFun,sqrt(tolf)};
  done = true;
elseif Delta < 2*eps
  EXITFLAG = -3;
  msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve', ...
      Delta,'',2*eps};
  done = true;
elseif iter >= maxit
  EXITFLAG = 0;
  msgData = {'trustnleqn',10,verbosity > 0,detailedExitMsg,'fsolve', ...
      [],optionFeedback.MaxIter,maxit};
  done = true;
elseif numFevals >= maxfunc
  EXITFLAG = 0;
  msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve', ...
      [],optionFeedback.MaxFunEvals,maxfunc};
  done = true;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Delta = updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                             alpha1,alpha2,DeltaMax,faultTolStruct)
%updateDelta updates the trust region radius in TRUSTNLEQN.
%
%   updateDelta updates the trust region radius based on the value of
%   ratio and the norm of the scaled step.

if ~faultTolStruct.currTrialWellDefined
    % Shrink the trust region if any element of the function vector
    % at the new step is not defined (i.e. it's inf/NaN/complex). The
    % update for Delta in this case follows that in snls.
    Delta = min(normdscal/20,Delta/20);
elseif ratio < eta1
    Delta = alpha2*normdscal;
elseif ratio >= eta2
    Delta = max(Delta,alpha1*normdscal);
end
Delta = min(Delta,DeltaMax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scalMat = getscalMat(nvar,JAC,scalemin,scalemax)
%getscalMat computes the scaling matrix in TRUSTNLEQN.
%
%   getscalMat computes the scaling matrix based on the norms 
%   of the columns of the Jacobian.

scalMat = ones(nvar,1);
for i=1:nvar
  scalMat(i,1) = norm(JAC(:,i));
end
scalMat(scalMat<scalemin) = scalemin;  % replace small entries
scalMat(scalMat>scalemax) = scalemax;  % replace large entries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state,iter,numFevals, ...
    Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter','interrupt', or 'done'. 
%
% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = numFevals;
optimValues.fval = Fvec;
optimValues.stepsize = normd; 
optimValues.gradient = grad; 
optimValues.firstorderopt = normgradinf;
optimValues.trustregionradius = Delta;
optimValues.stepaccept = stepAccept;

xOutputfcn(:) = xvec;  % Set xvec to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init','interrupt'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('optim:trustnleqn:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
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
            error(message('optim:trustnleqn:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end
%--------------------------------------------------------------------------
function [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT sets the outputs arguments to be the values at the last call
% of the outputfcn during an 'iter' call (when these values were last known to
% be consistent). 

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn; 
Fvec = optimValues.fval;
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'trust-region-dogleg';
OUTPUT.firstorderopt = optimValues.firstorderopt; 
JAC = []; % May be in an inconsistent state



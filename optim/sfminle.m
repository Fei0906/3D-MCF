function[x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminle(funfcn,x,A,b,verb,options,defaultopt,...
    computeLambda,initialf,initialGRAD,initialHESS,Hstr,detailedExitMsg,computeConstrViolForOutput,varargin)
%

%SFMINLE Nonlinear minimization with linear equalities.
%
% Locate a local minimizer to 
%
%               min { f(x) :  Ax = b }.
%
%	where f(x) maps n-vectors to scalars.
%

%   Copyright 1990-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.22 $  $Date: 2011/08/29 20:34:13 $

caller = funfcn{2};
% Set algorithm in order to pupulate output.algorithm
algorithmName = 'trust-region-reflective';

% Check to see if the function value at the initial point is well-defined.
% If not, then terminate immediately.
if ~isfinite(initialf) || ~isreal(initialf) 
    error(message('optim:sfminle:UsrObjUndefAtX0', caller));    
end

%   Initialization
xcurr = x(:); % x has "the" shape; xcurr is a vector
numFunEvals = 1;  % done in calling function fmincon

nVar = length(xcurr); 
iter = 0;
header = sprintf(['\n                                Norm of      First-order \n',...
        ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstrFirstIter = ' %5.0f      %13.6g                  %12.3g                ';
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';

if nVar == 0
    error(message('optim:sfminle:InvalidN'))
end

% Add Preconditioner to defaultopt. It's not a documented option, and is
% not added to defaultopt at the user-facing function level.
defaultopt.Preconditioner = @preaug;

% Read in user options
pcmtx = optimget(options,'Preconditioner',defaultopt,'fast');
mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
% Check if name clash
functionNameClashCheck('HessMult',mtxmpy,'hessMult_optimInternal', ...
    'optim:sfminle:HessMultNameClash');

% Use internal Hessian-multiply function if user does not provide HessMult function 
% or options.Hessian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungradhess') && ~strcmpi(funfcn{1},'fun_then_grad_then_hess'))
    mtxmpy = @hessMult_optimInternal;
end

pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tol2 = optimget(options,'TolX',defaultopt,'fast') ;
tol1 = optimget(options,'TolFun',defaultopt,'fast') ;
maxiter = optimget(options,'MaxIter',defaultopt,'fast') ;
maxfunevals = optimget(options,'MaxFunEvals',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast') ;  % pcgtol = .1;
pcgtol = min(pcgtol, 1e-2);  % better default
kmax = optimget(options,'MaxPCGIter',defaultopt,'fast') ;
if ischar(kmax)
    if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
        kmax = max(1,floor(nVar/2));
    else
        error(message('optim:sfminle:InvalidMaxPCGIter'))
    end
end

if ischar(maxfunevals)
    if isequal(lower(maxfunevals),'100*numberofvariables')
        maxfunevals = 100*nVar;
    else
        error(message('optim:sfminle:InvalidMaxFunEvals'))
    end
end

% Prepare strings to give feedback to users on options they have or have not set.
% These are used in the exit messages.
optionFeedback = createOptionFeedback(options);

% Handle the output
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot functions
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

dnewt = []; 
snod = [];
done = false; 
posdef = 1; npcg = 0;
pcgit = 0; delta = 100;nrmsx = 1; ratio = 0; degen = inf; 
dv = ones(nVar,1);
oval = inf;  gradf = zeros(nVar,1); newgrad = gradf; Z = []; 

%   Remove (numerical) linear dependencies from A
AA = A; bb = b;
[A,b] = dep(AA,[],bb);
m = size(A,1); mm = size(AA,1);
if verb > 2 && m ~= mm
    fprintf(getString(message('optim:sfminle:RemovedLinConstrDepend',mm-m)));
end

%   Get feasible: nearest feas. pt. to xstart
xcurr = feasibl(A,b,xcurr);

% Make xcurr conform to the user's input x
x(:) = xcurr;

if ~isempty(Hstr) % use sparse finite differencing
    
    switch funfcn{1}
        case 'fungrad'
            val = initialf; gradf(:) = initialGRAD;
        case 'fun_then_grad'
            val = initialf; gradf(:) = initialGRAD;
        otherwise
            error(message('optim:sfminle:UndefinedCalltypeInFMINCON'))
    end
    
    %      Determine coloring/grouping for sparse finite-differencing
    p = colamd(Hstr)'; p = (nVar+1)*ones(nVar,1)-p; group = color(Hstr,p);
    H = sfd(x,gradf,Hstr,group,[],options.DiffMinChange,options.DiffMaxChange, ...
            funfcn,varargin{:});
    
else % user-supplied computation of H or dnewt
    switch funfcn{1}
        case 'fungradhess'
            val = initialf; gradf(:) = initialGRAD; H = initialHESS;
        case 'fun_then_grad_then_hess'
            val = initialf; gradf(:) = initialGRAD; H = initialHESS;
        otherwise
            error(message('optim:sfminle:UndefinedCalltypeInFMINCON'))
    end
end
[~,pnewt] = size(gradf);

%   Extract the Newton direction?
if pnewt == 2, dnewt = gradf(1:nVar,2); end
PT = findp(A);
[g,LZ,UZ,pcolZ,PZ] = project(A,-gradf(1:nVar,1),PT);
gnrm = norm(g,inf);

if verb > 2
    disp(header)
end

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'init',iter,numFunEvals, ...
        val,[],[],[],pcgit,[],[],[],delta,A,b,varargin{:});
    if stop
        [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
            cleanUpInterrupt(xOutputfcn,optimValues,npcg,verb,detailedExitMsg,algorithmName,caller);
        return;
    end
end

% Create the fault tolerance structure
faultTolStruct = createFaultTolStruct(true);
newFaultTolStruct = faultTolStruct;

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(iter) S.T. f(xcurr(iter)) IS DECREASING.
while ~done
    if any(~isfinite(gradf))
        error(message('optim:sfminle:InvalidUserFunction', caller))
    end
    
    % Display
    if verb > 2
        if iter > 0
            if faultTolStruct.undefObj 
                fprintf(getString(message('optimlib:commonMsgs:ObjInfNaNComplex', ...
                    faultTolStruct.undefValue)));
            end                        
            currOutput = sprintf(formatstr,iter,val,nrmsx,gnrm,pcgit);            
        else
            currOutput = sprintf(formatstrFirstIter,iter,val,gnrm);
        end
        disp(currOutput);
    end
    
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'iter',iter,numFunEvals, ...
            val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,A,b,varargin{:});
        if stop  % Stop per user request.
            [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
                cleanUpInterrupt(xOutputfcn,optimValues,npcg,verb,detailedExitMsg,algorithmName,caller);
            return;
        end
    end
    
    %     TEST FOR CONVERGENCE
    diff = abs(oval-val); 
    oval = val; 
    if gnrm < tol1 && posdef == 1
        done = true; 
        EXITFLAG = 1;
        if iter == 0
            msgFlag = 100;
        else
            msgFlag = EXITFLAG;
        end
        % Call createExitMsg with createExitMsgExitflag = 100 if x0 is
        % optimal, otherwise createExitMsgExitflag = 1
        outMessage = createExitMsg('sfminle',msgFlag,verb > 1,detailedExitMsg,caller, ...
            gnrm,optionFeedback.TolFun,tol1);

    elseif (nrmsx < .9*delta) && (ratio > .25) && (diff < tol1*(1+abs(oval)))
        done = true; 
        EXITFLAG = 3;
        outMessage = createExitMsg('sfminle',EXITFLAG,verb > 1,detailedExitMsg,caller, ...
            diff/(1+abs(oval)),optionFeedback.TolFun,tol1);        
    elseif iter > 1 && nrmsx < tol2
        done = true; 
        EXITFLAG = 2;
        outMessage = createExitMsg('sfminle',EXITFLAG,verb > 1,detailedExitMsg,caller, ...
            nrmsx,optionFeedback.TolX,tol2);
    elseif iter > maxfunevals
        done = true; 
        EXITFLAG = 0;
        outMessage = createExitMsg('sfminle',EXITFLAG,verb > 0,detailedExitMsg,caller, ...
            [],optionFeedback.MaxFunEvals,maxfunevals);
    elseif iter > maxiter
        done = true; 
        EXITFLAG = 0;
        % Call createExitMsg with createExitMsgExitflag = 10 for MaxIter exceeded
        outMessage = createExitMsg('sfminle',10,verb > 0,detailedExitMsg,caller, ...
            [],optionFeedback.MaxIter,maxiter);
    end
    
    % Reset fault tolerance structure
    faultTolStruct = newFaultTolStruct;
 
    %     Step computation
    if ~done
        % OutputFcn call
        if haveoutputfcn % Call output functions (we don't call plot functions with 'interrupt' flag)
            [~, ~, stop] = callOutputAndPlotFcns(outputfcn,{},xcurr,xOutputfcn,'interrupt',iter,numFunEvals, ...
                val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,A,b,varargin{:});
            if stop  % Stop per user request.
                [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
                    cleanUpInterrupt(xOutputfcn,optimValues,npcg,verb,detailedExitMsg,algorithmName,caller);
                return;
            end
        end

        %       Determine trust region correction
        sx = zeros(nVar,1); 
        oposdef = posdef;
        [sx,snod,qp,posdef,pcgit,Z] = trdg(xcurr,gradf(:,1),H,...
            delta,g,mtxmpy,pcmtx,pcflags,...
            pcgtol,kmax,A,zeros(m,1),Z,dnewt,options,defaultopt,...
            PT,LZ,UZ,pcolZ,PZ,varargin{:});
        
        if isempty(posdef), 
            posdef = oposdef; 
        end
        nrmsx=norm(snod); 
        npcg=npcg + pcgit;
        newx=xcurr + sx; 
        
        % Make newx conform to user's input x
        x(:) = newx;
        %       Evaluate f,g,  and H
        if ~isempty(Hstr) % use sparse finite differencing
            switch funfcn{1}
                case 'fungrad'
                    [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
                case 'fun_then_grad'
                    newval = feval(funfcn{3},x,varargin{:}); 
                    newgrad(:) = feval(funfcn{4},x,varargin{:});
                otherwise
                    error(message('optim:sfminle:UndefinedCalltypeInFMINCON'))
            end
            
            newH = sfd(x,newgrad,Hstr,group,[],options.DiffMinChange,options.DiffMaxChange, ...
                       funfcn,varargin{:});
            
        else % user-supplied computation of H or dnewt
            switch funfcn{1}
                case 'fungradhess'
                    [newval,newgrad(:),newH] = feval(funfcn{3},x,varargin{:});
                case 'fun_then_grad_then_hess'
                    newval = feval(funfcn{3},x,varargin{:}); 
                    newgrad(:) = feval(funfcn{4},x,varargin{:});
                    newH = feval(funfcn{5},x,varargin{:});
                otherwise
                    error(message('optim:sfminle:UndefinedCalltypeInFMINCON'))
            end
            
        end
        numFunEvals = numFunEvals + 1;

        % Update fault tolerance structure
        faultTolStruct = updateFaultTolStruct(faultTolStruct,...
            newval, verb > 2);
        
        % Update trust region radius
        if ~faultTolStruct.currTrialWellDefined
            % Shrink the trust region if the function value at the new step
            % is not defined (i.e. it's inf/NaN/complex)
            delta = min(nrmsx/20,delta/20);
        else
            % Update trust region using change in objective and in
            % quadratic model            
            [~,pnewt] = size(newgrad);
            if pnewt == 2,
                dnewt = newgrad(1:nVar,2);
            end
            aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
            ratio = (newval + aug -val)/qp;
            if (ratio >= .75) && (nrmsx >= .9*delta),
                delta = 2*delta;
            elseif ratio <= .25,
                delta = min(nrmsx/4,delta/4);
            end
        end
        
        % Accept the step if the function value is defined at the new step
        % and the function value is lower than at the current step
        if faultTolStruct.currTrialWellDefined && newval < val
            xcurr = newx; 
            val = newval; 
            gradf = newgrad; 
            H = newH;
            Z = [];
            if pnewt == 2, 
                dnewt = newgrad(1:nVar,2); 
            end
            g = project(A,-gradf(:,1),PT,LZ,UZ,pcolZ,PZ);
            gnrm = norm(g,inf);
            
            %          Extract the Newton direction?
            if pnewt == 2, 
                dnewt = newgrad(1:nVar,2); 
            end
        end
        iter = iter+1;
        
    end % if ~done
end % while

if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'done',iter,numFunEvals, ...
        val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,A,b,varargin{:});
end


HESSIAN = H;
GRAD = g;
FVAL = val;
if computeLambda
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
    warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
    warningstate2 = warning('off', 'MATLAB:singularMatrix');

    LAMBDA.eqlin = -A'\gradf;

    % Restore the warning states to their original settings
    warning(warningstate1)
    warning(warningstate2)
    LAMBDA.eqnonlin = zeros(0,1); % An empty column-vector
    LAMBDA.ineqlin = zeros(0,1);  % An empty column-vector
    % Assign zero multipliers for +/- Inf upper and lower bounds
    % This is consistent with other algorithms in fmincon
    LAMBDA.lower = zeros(nVar,1);
    LAMBDA.upper = zeros(nVar,1);
    LAMBDA.ineqnonlin = zeros(0,1); % An empty column-vector
else
    LAMBDA = [];
end
OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = algorithmName;
OUTPUT.message = outMessage;
if computeConstrViolForOutput
    OUTPUT.constrviolation = norm(A*xcurr-b, inf);
else
    OUTPUT.constrviolation = [];
end
x(:) = xcurr;
%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state,iter,numFunEvals, ...
    val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,A,b,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then
% calls the outputfcn/plotfcns.  
%
% state - can have the values 'init','iter','interrupt', or 'done'. 

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = numFunEvals;
optimValues.fval = val;
optimValues.stepsize = nrmsx;
optimValues.gradient = g;
optimValues.firstorderopt = gnrm;
optimValues.cgiterations = pcgit; 
optimValues.positivedefinite = posdef;
optimValues.ratio = ratio;
optimValues.degenerate = min(degen,1);
optimValues.trustregionradius = delta;
optimValues.procedure = '';
optimValues.constrviolation = norm(A*xvec-b, inf);

xOutputfcn(:) = xvec;  % Set x to have user expected size
stop = false;
if ~isempty(outputfcn)
    switch state
        case {'iter','init','interrupt'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('optim:sfminle:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
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
            error(message('optim:sfminle:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end

%--------------------------------------------------------------------------
function [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = cleanUpInterrupt(xOutputfcn,optimValues,npcg, ...
                                                    verb,detailedExitMsg,algorithmName,caller)
% CLEANUPINTERRUPT updates or sets all the output arguments of SFMINBX when the optimization 
% is interrupted.  The HESSIAN and LAMBDA are set to [] as they may be in a
% state that is inconsistent with the other values since we are
% interrupting mid-iteration.

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn;
FVAL = optimValues.fval;
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.stepsize = optimValues.stepsize;
OUTPUT.algorithm = algorithmName; 
OUTPUT.firstorderopt = optimValues.firstorderopt; 
OUTPUT.cgiterations = npcg; % total number of CG iterations
OUTPUT.message = createExitMsg('sfminle',EXITFLAG,verb > 0,detailedExitMsg,caller);
OUTPUT.constrviolation = optimValues.constrviolation;
GRAD = optimValues.gradient;
HESSIAN = []; % May be in an inconsistent state
LAMBDA = []; % May be in an inconsistent state


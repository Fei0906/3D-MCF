function [options,optionFeedback] = getIpOptions(options,nVar,mEq,nonlconflag,defaultopt,defaultHessMemory,defaultPivotThreshold)
%

%getIpOptions Helper function that reads options based on the values in
% user-supplied structure "options" and on defaultopt

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/05/09 01:05:54 $

% Prepare strings to give feedback to users on options they have or have not set.
% These are used in the exit messages.
optionFeedback = createOptionFeedback(options);

% Read in all options fields that don't need processing
fieldsThatDontNeedProcessing = {'AlwaysHonorConstraints'; ...
    'Diagnostics';'Display'; ...
    'FunValCheck';'HessFcn';'HessMult';'InitBarrierParam'; ...
    'MaxIter';'ObjectiveLimit';'OutputFcn';'PlotFcns';'ScaleProblem'; ...
    'TolCon';'TolFun';'TolGradCon';'TolProjCG';'TolProjCGAbs';'TolX';'UseParallel'};

for i = 1:numel(fieldsThatDontNeedProcessing)
    field = fieldsThatDontNeedProcessing{i};
    options.(field) = optimget(options,field,defaultopt,'fast');
end

% Now read-in all fields that can be either a string or a numeric value
% Each row of cell array below contains: 
% field name, allowed string, equivalent numeric value, datatype (with appropriate article)
% for error message
numericOrStringFields = ... % NB: This is a *2-dimensional* cell-array
    {'MaxFunEvals','100*numberofvariables',100*nVar,'an integer';
    'InitTrustRegionRadius','sqrt(numberofvariables)',sqrt(nVar),'a real';
    'MaxProjCGIter','2*(numberofvariables-numberofequalities)',max(2*(nVar-mEq),0),'a matrix'};

for i = 1:size(numericOrStringFields,1)
    field = numericOrStringFields{i,1};
    allowedString = numericOrStringFields{i,2};
    equivNumValue = numericOrStringFields{i,3};
    dataType = numericOrStringFields{i,4};
    [options.(field),ME] = getNumericOrStringFieldValue(field,allowedString, ...
        equivNumValue,dataType,options,defaultopt);
    if ~isempty(ME)
        throwAsCaller(ME)
    end
end

% Now read-in fields that need special processing

% Option Hessian can be either a string or a cell array. So as not
% to have to check for this inside barrier.m at every iteration,
% we create one or two internal options: HessType and HessMemory.
Hessian = optimget(options,'Hessian',defaultopt,'fast');
options.HessMemory = defaultHessMemory;
if ischar(Hessian)
    % If Hessian user-supplied, on, fin-diff-grads, gradients must be on
    if ( strcmpi(Hessian,'user-supplied') || strcmpi(Hessian,'on') ...
            || strcmpi(Hessian,'fin-diff-grads') ) && ...
            ( strcmpi(options.GradObj,'off') || nonlconflag && strcmpi(options.GradConstr,'off') )
        mexcptn = MException('optim:getIpOptions:HessButNoGrads', ...
            'Must set GradObj = ''on'' and (if nonlinear constraints are present)\nGradConstr = ''on'' in order to use option Hessian = ''%s''.',Hessian);
        throwAsCaller(mexcptn);
    end    
    if strcmpi(Hessian,'off')
        % No Hessian function provided and old value; set to default
        options.HessType = defaultopt.Hessian;
    elseif strcmpi(Hessian,'user-supplied') || strcmpi(Hessian,'on')
        % Hessian provided, either via options HessFcn or HessMult
        if ~isempty(options.HessFcn) 
            options.HessType = 'user-supplied';
            if ~isempty(options.HessMult)
                % Both HessFcn and HessMult provided; honor HessFcn
                warning(message('optim:getIpOptions:BothHessFcnHessMult'))
            end
        elseif ~isempty(options.HessMult)
            options.HessType = 'hessmult';
        else 
            % Hessian user-supplied but no HessFcn nor HessMult
            mexcptn = MException('optim:getIpOptions:BadHessOptions', ...
                'Hessian option set to ''%s'' but no Hessian function provided in options HessFcn nor in HessMult.',Hessian);
            throwAsCaller(mexcptn);
        end
    else % any of the current legal string values for Hessian option
        options.HessType = Hessian;
    end
elseif iscell(Hessian)
    options.HessType = Hessian{1};
    if numel(Hessian) > 1
        options.HessMemory = Hessian{2};
    end
end

% SubproblemAlgorithm corresponds to three internal options: IpAlgorithm,
% LinearSystemSolver, and PivotThreshold.
% The value of SubproblemAlgorithm can be either a string or a cell array
SubproblemAlgorithm = optimget(options,'SubproblemAlgorithm',defaultopt,'fast');

% Store subproblem algorithm in a string for easy processing
% and read-in pivot threshold
% Because of error checking in options functions, SubproblemAlgorithm
% can only be either a string or a cell array.
if ischar(SubproblemAlgorithm)
    SubproblemAlgorithm_string = SubproblemAlgorithm;
    options.PivotThreshold = defaultPivotThreshold; 
elseif iscell(SubproblemAlgorithm)
    SubproblemAlgorithm_string = SubproblemAlgorithm{1};
    if numel(SubproblemAlgorithm) > 1
        % Pivot threshold passed in with option
        options.PivotThreshold = SubproblemAlgorithm{2};
    else
        options.PivotThreshold = defaultPivotThreshold; 
    end
end

if strcmpi(SubproblemAlgorithm_string,'ldl-factorization')
    options.IpAlgorithm = 'direct';
    options.LinearSystemSolver = 'ldl-factorization';
elseif strcmpi(SubproblemAlgorithm_string,'cg') 
    options.IpAlgorithm = 'cg';
    options.LinearSystemSolver = 'ldl-factorization';
end

% Must select SubproblemAlgorithm = 'cg' in order to user fin-diff-grads or HessMult
if ~strcmpi(options.IpAlgorithm,'cg')
    if strcmpi(options.HessType,'fin-diff-grads')
        mexcptn = MException('optim:getIpOptions:FinDiffGradsAndNotCG', ...
            'In order to use option Hessian = ''fin-diff-grads'', option SubproblemAlgorithm\n must be set to ''cg''.');
        throwAsCaller(mexcptn);
    elseif strcmpi(options.HessType,'hessmult')
        mexcptn = MException('optim:getIpOptions:HessMultAndNotCG', ...
            'In order to use option HessMult, option SubproblemAlgorithm must be set to ''cg''.');
        throwAsCaller(mexcptn);
    end
end

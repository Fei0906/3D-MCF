function resultsStruct = createResultsStruct(solverName,useValues)
%CREATERESULTSSTRUCT Create results structure for different solvers
%   Create result structure for 'solverName'. The optional second argument 
%   is used to populate the result structure 'resultStruct' with the values 
%   from 'useValues'.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2010/10/08 17:13:37 $

if nargin < 2
    useValues = [];
end
% The fields in the structure are in the same order as they are returned by
% the corresponding solver. 
switch solverName
    case 'fmincon' %1
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[],'lambda',[],'grad',[],'hessian',[]);
    case 'fminunc' %2
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[],'grad',[],'hessian',[]);
    case 'lsqnonlin' %3
        resultsStruct = struct('x',[],'resnorm',[],'residual',[], ...
            'exitflag',[],'output',[],'lambda',[],'jacobian',[]);
    case 'lsqcurvefit' %4
        resultsStruct = struct('x',[],'resnorm',[],'residual',[], ...
            'exitflag',[],'output',[],'lambda',[],'jacobian',[]);
    case 'fsolve' %5
    resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[],'jacobian',[]);
    case 'lsqlin' %6
            resultsStruct = struct('x',[],'resnorm',[],'residual',[], ...
                'exitflag',[],'output',[],'lambda',[]);
    case 'fgoalattain' %7
            resultsStruct = struct('x',[],'fval',[],'attainfactor',[], ...
                'exitflag',[],'output',[],'lambda',[]);
    case 'fminimax' %8
            resultsStruct = struct('x',[],'fval',[],'maxfval',[], ...
                'exitflag',[],'output',[],'lambda',[]);   
    case 'fseminf' %9
            resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
                'output',[],'lambda',[]);
    case 'linprog' %10
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[],'lambda',[]);
    case 'quadprog' %11
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[],'lambda',[]);
    case 'bintprog' %12
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[]);
    case 'fminbnd' %13
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[]);
    case 'fminsearch' %14
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[]);
    case 'fzero' %15
        resultsStruct = struct('x',[],'fval',[],'exitflag',[], ...
            'output',[]);
    case 'lsqnonneg' %16
        resultsStruct = struct('x',[],'resnorm',[],'residual',[], ...
            'exitflag',[],'output',[],'lambda',[]);
    case 'ga' %17
        resultsStruct = struct('x',[],'fval',[],'exitflag',[],'output',[], ...
            'population',[],'score',[]);
    case 'patternsearch' %18
        resultsStruct = struct('x',[],'fval',[],'exitflag',[],'output',[]);
    case 'simulannealbnd' %19
        resultsStruct = struct('x',[],'fval',[],'exitflag',[],'output',[]);
    case 'gamultiobj' %20
        resultsStruct = struct('x',[],'fval',[],'exitflag',[],'output',[], ...
          'population',[],'score',[]);
    otherwise
        error(message('optim:createResultsStruct:UnrecognizedSolver'));
end

% Copy the values from the struct 'useValues' to 'resultsStruct'
if ~isempty(useValues)
    copyfields = fieldnames(resultsStruct);
    Index = ismember(fieldnames(resultsStruct),fieldnames(useValues))
    for i = 1:length(Index)
        if Index(i)
            resultsStruct.(copyfields(i)) = useValues.(copyfields(i));
        end
    end
end

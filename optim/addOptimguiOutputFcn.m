function optStruct = addOptimguiOutputFcn(optStruct,solverName)
%

%addOptimguiOutputFcn private to OPTIMTOOL.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/05/09 01:05:43 $

% Determine the appropriate field and value for the solver
switch solverName
    case {'fmincon','fminunc','lsqnonlin','lsqcurvefit','linprog', ...
            'quadprog','bintprog','fgoalattain','fminimax','fseminf', ...
            'fminsearch','fzero','fminbnd','fsolve','lsqlin','lsqnonneg'}
        outputFcnFieldName = 'OutputFcn';
        outputFcnValue = @optimtooloutput;
    case 'ga'
        outputFcnFieldName = 'OutputFcns';
        outputFcnValue = @gatooloutput;
    case 'gamultiobj'
        outputFcnFieldName = 'OutputFcns';
        outputFcnValue = @gamultiobjtooloutput;
    case 'patternsearch'
        outputFcnFieldName = 'OutputFcns';
        outputFcnValue = @psearchtooloutput;
    case 'simulannealbnd'
        outputFcnFieldName = 'OutputFcns';
        outputFcnValue = @satooloutput;
end

% We need to adjust 'OutputFcn(s)' field so that the output function is
% called in every iteration. The GUI and solvers interact through the output 
% function.
if isempty(optStruct.(outputFcnFieldName))
    optStruct.(outputFcnFieldName) = outputFcnValue;
elseif iscell(optStruct.(outputFcnFieldName))
    % Add 'outputFcnValue' at the end of the array
    optStruct.(outputFcnFieldName){end+1} = outputFcnValue;
else % Make it a cell
    optStruct.(outputFcnFieldName) = {optStruct.(outputFcnFieldName)};
    % Add 'outputFcnValue' as output function
    optStruct.(outputFcnFieldName){end+1} = outputFcnValue;
end

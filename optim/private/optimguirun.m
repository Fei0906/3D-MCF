function [err,x,fval,exitMessage,nrow,ncolx,ncolf] = optimguirun(hashProb, hashOpt)
%OPTIMGUIRUN Optimization Toolbox GUI 'Start' button callback function.
%   Two arguments 'hashProb' and 'hashOpt' are Java hash tables for 
%   problem model and options model respectively.
%   The output 'err' is the error string returned by either readOptimHashTable
%   or callSolver functions. x,fval,exitMessage are outputs from the solver
%   and [nrow,ncol] is the size of 'x' which is needed to display 'x' in the GUI. 

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2011/11/09 16:56:08 $

err = '';
x = '';
fval = '';
exitMessage = '';
% Size of the result vector 'X' and 'fval' (nrows is always the number of
% solution)
nrow = [];
ncolx = [];
ncolf = [];
% Maximum length of 'X' vector to be shown in the GUI
MAX_NUM_ELEMENT_SHOW = 100;
% Get modified fields from the GUI
[probStruct,optStruct,errProb,errOpt] = readOptimHashTable(hashProb,hashOpt);
if ~isempty(errProb)
    err = getString(message('optim:optimtool:RuntimeErrorHeaderProblem',errProb));
    return;
elseif ~isempty(errOpt)
    err = getString(message('optim:optimtool:RuntimeErrorHeaderOptions',errOpt));
    return;
end
% Add appropriate output function for the 'solver' so
% that the GUI and the solver can interact in each iteration
optStruct = addOptimguiOutputFcn(optStruct,probStruct.solver);
% Store current warning state
[lastmsg, lastid] = lastwarn;
lastwarn('');

try
    % Call solver and save the result structure to MATLAB workspace (appdata)
    resultStruct = callSolver(probStruct,optStruct);
    setappdata(0,'optimTool_results_121677',resultStruct);
    % Retrieve the basic parts of the exit message
    [summary,defaultBody] = optimExitMsgParts(resultStruct.output.message);
    % Form simple message from retrieved parts
    exitMessage = sprintf('%s\n\n%s\n',summary,defaultBody);
    x = resultStruct.x;
    
    % Set iteration number in the GUI
    optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
    if ~isempty(optimtoolGui)
        % Update iteration number in the GUI
        try
            javaMethodEDT('setIteration',optimtoolGui,value2RHS(resultStruct.output.iterations));
        catch ME % GA solvers have 'generations' and not 'iterations'
            javaMethodEDT('setIteration',optimtoolGui,value2RHS(resultStruct.output.generations));
        end
    end

    % Least square solvers have 'resnorm' and not 'fval'
    try
        fval = resultStruct.fval;
    catch ME
        fval = resultStruct.resnorm;
    end
    % Return argument fval
    if ndims(fval) < 3 && numel(fval) <= MAX_NUM_ELEMENT_SHOW
        if isscalar(fval) % Single objective
            ncolf = 0; % GUI expects one less
        else % Multiple solution (only gamultiobj)
          [~,ncolf] = size(fval);
        end
    else
        ncolf = -1;
        fval  = [];
    end
    % Return argument x
    if ndims(x) < 3 && (isnumeric(x) || isa(x,'double')) && ...
            numel(x) <= MAX_NUM_ELEMENT_SHOW
            [nrow, ncolx] = size(x);          
    else
        nrow  = -1;
        ncolx = -1;
        x = [];
    end
    % Save the random states if it is returned in the resultStruct
    if isfield(resultStruct.output,'rngstate')
        probStruct.rngstate = resultStruct.output.rngstate;
        setoptimrandstates(probStruct,java.util.Hashtable,false);
    end
catch ME
    nrow = -1;
    ncolx = -1;
    ncolf = -1;
    x = [];
    err = ME.message;
    % Remove html tags from the error message
    err = regexprep(err,'</?(\w+).*?>','');
    % Set 'Error running optimization' message in the GUI
    endMessage = getString(message('optim:optimtool:OptimtoolStatusRuntimeErr'));
    optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
    if ~isempty(optimtoolGui)
        javaMethodEDT('appendResults',optimtoolGui,endMessage);
    end    
end

% Restore the warning if there were no warnings from the solver
if isempty(lastwarn)
    lastwarn(lastmsg,lastid);
end


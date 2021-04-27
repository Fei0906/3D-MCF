function results = callSolver(probStruct,optStruct)
%callSolver Call appropriate solver and return the results structure.
%   The input arguments 'probStruct' and 'optStruct' are problem and options structure.
%   The output 'results' is a structure containing all the outputs from solver.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2011/06/30 16:37:02 $

solver = probStruct.solver;
% Update the problem structure 'probStruct' for the solver
probStruct = createProblemStruct(solver,[],probStruct);
% Add options to the problem structure (which is a required field)
probStruct.options = createOptionsStruct(solver,optStruct); 
% Create result structure for the solver
results = createResultsStruct(solver);
resultsFields = fieldnames(results);
numfields = length(resultsFields);
% Initialize cell array to capture the output from the solver
solverOutput = cell(1,numfields);
% Set 'optimization running' message in the GUI
startMessage = getString(message('optim:optimtool:OptimtoolStatusRunning'));
optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
if ~isempty(optimtoolGui)
    javaMethodEDT('appendResults',optimtoolGui,startMessage);
end
% Call solver and put results in the structure
[solverOutput{:}] = feval(str2func(solver),probStruct);
for i = 1:numfields
    results.(resultsFields{i}) = solverOutput{i};
end

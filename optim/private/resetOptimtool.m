function [selection, optionsModel] = resetOptimtool()
% RESETOPTIMTOOL Reset OPTIMTOOL
%   Private to OPTIMTOOL

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:44 $

% createProblemStruct first two arguments are solverName and defaultSolver.
% In this case, we want to set solverName to "all" so that all fields will
% be cleared. The "defaultSolver" is fmincon

probStruct = createProblemStruct('all', 'fmincon'); 
setappdata(0,'optimTool_Problem_Data',probStruct);

% Create the empty options structure
options = createOptionsStruct('all',struct('Display', 'off')); % This is the default for the GUI
setappdata(0,'optimTool_Options_Data',options);

% Remove appdata structures used by optimtool
if isappdata(0,'optimTool_results_121677')
    rmappdata(0,'optimTool_results_121677');
end

% Reset hash table for problem and options
if isappdata(0,'optimTool_Problem_HashTable')
    rmappdata(0,'optimTool_Problem_HashTable');
end

if isappdata(0,'optimTool_Options_HashTable')
    rmappdata(0,'optimTool_Options_HashTable');
end

[selection, optionsModel] = optimguiImportOptions(1);

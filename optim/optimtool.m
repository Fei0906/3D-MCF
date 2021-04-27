function optimtool(inputArg)
%OPTIMTOOL Optimization Toolbox Graphical User Interface.
%   OPTIMTOOL starts a graphical user interface of the optimization
%   solvers from MATLAB, Optimization Toolbox, and if a license is
%   available, Global Optimization Toolbox. OPTIMTOOL can be used for 
%   editing the default options, selecting, and running solver.
%
%   OPTIMTOOL(INPUTARG) starts the Optimization Tool with INPUTARG.
%   INPUTARG can be either optimization options structure or optimization
%   problem structure. An options structure can be created using OPTIMSET
%   or by using the export option from OPTIMTOOL. See OPTIMSET for a
%   detailed description on creating options. A problem structure can be
%   created in OPTIMTOOL and exported to MATLAB workspace. INPUTARG can
%   also be a name of optimization solver (character string). OPTIMTOOL
%   will open with the default options and problem fields for the solver
%   INPUTARG.
%
%   You can import optimization options structure to OPTIMTOOL from MATLAB
%   workspace and modify it in OPTIMTOOL. You can also export an options
%   structure to MATLAB workspace from OPTIMTOOL. You also can import a
%   problem structure from MATLAB workspace to OPTIMTOOL. You can export a
%   problem structure to MATLAB workspace from OPTIMTOOL. See help for
%   optimization solver for a detailed description of specific fields in a
%   problem structure.
%
%   OPTIMTOOL can be used to run all optimization solvers after setting a
%   problem. You can export the results e.g., X and FVAL, etc as a result
%   structure to the MATLAB workspace.
%
%   See also OPTIMSET.

%   Copyright 2006-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.12 $  $Date: 2011/06/30 16:36:56 $

% Can't proceed unless we have java support
if ~usejava('swing')
    error(message('optim:optimtool:missingJavaSwing'));
end
% Default solver is fmincon
defaultSolver = 'fmincon';
% Required field names for problem structure
requiredFields = {'solver','options'};
validValues = {fieldnames(createProblemStruct('solvers')), {} };
probFieldnames = fieldnames(createProblemStruct('all',defaultSolver));
% Fieldnames of options structure
optionsFieldnames = fieldnames(createOptionsStruct('all'));
% We check for one matching field
minNumberOfOptionsToCheck = 1; % Display
optimGUI = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI');
if nargin < 1 % No input arguments
    if isempty(optimGUI)
        % Create the empty problem structure and the options structure for
        % the default solver. Also, create Java hash table for these structures.
        options = createOptionsStruct('all',struct('Display', 'off')); % This is the default for the GUI
        optionsModel = createHashTable(options);
        probStruct = createProblemStruct('all',defaultSolver);
        problemModel = createHashTable(probStruct);
    else
       javaMethodEDT('toFront',optimGUI);
       return;
    end
elseif isstruct(inputArg) % Input argument is either a problem structure or options structure
    % Problem structure must have all the required fields
    [validProblemStruct,errmsg,inputArg] = validOptimProblemStruct(inputArg,requiredFields,validValues);
    if validProblemStruct
        % Make sure that the GUI is not open; else values will be replaced
        if ~isempty(optimGUI)
            button = questdlg(getString(message('optim:optimtool:QuestDlgReplaceCurrProbAndOpts')), ...
                getString(message('optim:optimtool:ReplaceCurrVals')), ... % Message
                getString(message('optim:optimtool:BtnYes')), ... % Yes
                getString(message('optim:optimtool:BtnNo')), ... % No
                getString(message('optim:optimtool:BtnNo'))); % Default choice: No
            if ~strcmp(button, getString(message('optim:optimtool:BtnYes')))
                javaMethodEDT('toFront',optimGUI);
                return;
            end
        end
        % Extract options structure from the problem and use 'optionsFieldnames' to create
        % a Java hash table. Note that the user supplied options structure may have spurious
        % fields which will be ignored. The fields we are interested in are in 'optionsFieldnames'.
        options = createOptionsStruct('all',inputArg.options); 
        inputArg = rmfield(inputArg,'options');
        optionsModel = createHashTable(options,optionsFieldnames,sprintf('%s.options',inputname(1)));
        % Similarly, we create the problem hash table.
        probStruct = createProblemStruct('all',defaultSolver,inputArg);
        problemModel = createHashTable(inputArg,probFieldnames,inputname(1)); % Create Java hashtable
        if ~validOptions(options,optionsFieldnames,minNumberOfOptionsToCheck)  % Is valid options structure?
            if ~isempty(inputname(1))
                msg = getString(message('optim:optimtool:InvalidInputOpts',inputname(1)));
            else
                msg = getString(message('optim:optimtool:InvalidOpts'));
            end
            errordlg(msg,'Optimization Tool');
            return;
        end
    elseif validOptions(inputArg,optionsFieldnames,minNumberOfOptionsToCheck)  % Is valid options structure?
        % Make sure that the GUI is not open; else values will be replaced
        if ~isempty(optimGUI)
            button = questdlg(getString(message('optim:optimtool:QuestDlgReplaceCurrProbAndOpts')), ...
                getString(message('optim:optimtool:ReplaceCurrVals')), ... % Message
                getString(message('optim:optimtool:BtnYes')), ... % Yes
                getString(message('optim:optimtool:BtnNo')), ... % No
                getString(message('optim:optimtool:BtnNo'))); % Default choice: No
            if ~strcmp(button, getString(message('optim:optimtool:BtnYes')))
                javaMethodEDT('toFront',optimGUI);
                return;
            end
        end

        options = createOptionsStruct('all',inputArg);
        optionsModel = createHashTable(options,optionsFieldnames,inputname(1));
        probStruct = createProblemStruct('all',defaultSolver);
        problemModel = createHashTable(probStruct);
    else % Error
        errordlg(errmsg,'Optimization Tool');
        return;
    end
elseif ischar(inputArg)
    if strcmpi(inputArg, 'close')
        if ~isempty(optimGUI)
            javaMethodEDT('dispose',optimGUI);
        end
        return;
    else % Must be a solver name
        % Make sure that the GUI is not open; else values will be replaced
        if ~isempty(optimGUI)
            button = questdlg(getString(message('optim:optimtool:QuestDlgOptimtoolAlreadyOpen')), ...
                getString(message('optim:optimtool:ReplaceCurrVals')), ... % Message
                getString(message('optim:optimtool:BtnYes')), ... % Yes
                getString(message('optim:optimtool:BtnNo')), ... % No
                getString(message('optim:optimtool:BtnNo'))); % Default choice: No
            if ~strcmp(button, getString(message('optim:optimtool:BtnYes')))
                javaMethodEDT('toFront',optimGUI);
                return;
            end
        end
        options = createOptionsStruct('all',struct('Display', 'off')); % This is the default for the GUI
        optionsModel = createHashTable(options);
        probStruct = createProblemStruct(inputArg,defaultSolver); % Create MATLAB problem structure
        problemModel = createHashTable(probStruct); % Create Java problem hashtable
    end
else
    error(message('optim:optimtool:invalidInput'));
end

% Save problem and options to appdata
setappdata(0,'optimTool_Problem_Data',probStruct);
setappdata(0,'optimTool_Options_Data',options);
% Reset Java hashtable which stores the changes
resetOptimtoolHashTable('optimTool_Problem_HashTable');
resetOptimtoolHashTable('optimTool_Options_HashTable');

if ~isempty(ver('distcomp'))
    hasParallel = true;
else
    hasParallel = false;
end

if isempty(optimGUI)
    if isappdata(0,'optimTool_results_121677') % Unique string for appdata
        rmappdata(0,'optimTool_results_121677');
    end
    solversName = fieldnames(createProblemStruct('solvers'));
    javaObjectEDT('com.mathworks.toolbox.optim.OptimGUI', problemModel, optionsModel,(sprintf('%s ',solversName{:})),hasParallel);    
else
    javaMethodEDT('updateOptimGUI',optimGUI,problemModel,optionsModel);
end


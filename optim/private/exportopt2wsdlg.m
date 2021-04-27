function  exportopt2wsdlg(hashProb, hashOpt)
%exportopt2wsdlg exports variables from optimtool to workspace. 
%   Possible variables are 'problem', 'options', and 'results' structures. This 
%   function presents a dialog with choices for exporting any or all of the three 
%   above structures. It also creates variables in MATLAB workspace. The output 'err' 
%   will be a non-empty string (containing error message, if any) when there is any 
%   error in reading 'hashProb' and 'hashOpt'.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2011/10/15 01:57:42 $

% Title for the dialog
title = getString(message('optim:optimtool:ExportOptsDlgTitle'));
% Properties of the dialog (setting units to pixels to guard against
% someone setting defaultfigureunits to something else)
hDialog = dialog('Visible', 'off', 'Name', title, 'WindowStyle', 'normal', ...
    'Unit', 'pixels', 'Tag', 'ExportToWSDlg');
% Default variable names for the structures to be exported
defaultVariableNames = {'optimproblem'; ' '; 'options'; 'optimresults'};
variableNames = createVarNames(defaultVariableNames);
% Cancel button for the dialog
cancelStr = getString(message('optim:optimtool:BtnCancel'));
cancelButton = uicontrol(hDialog,'String', cancelStr, 'Unit', 'pixels', ...
                                 'Callback', {@CancelCallback, hDialog}, ...
                                 'Tag', 'CancelBtn');
% OK button for the dialog                          
okStr = getString(message('optim:optimtool:BtnOK'));
okButton = uicontrol(hDialog,'String', okStr, 'Unit', 'pixels', ...
                             'Tag','OkBtn');

% Labels for choices given in the dialog
checkboxLabels = {getString(message('optim:optimtool:CheckBoxExportProbAndOpts')); ...
                  getString(message('optim:optimtool:CheckBoxExportInfoToResume'));...
                  getString(message('optim:optimtool:CheckBoxExportOpts'));...
                  getString(message('optim:optimtool:CheckBoxExportResults'))};
        
% Retrieve problem and options structure from the Java hash table
[probStruct,optStruct] = readOptimHashTable(hashProb, hashOpt);

% Retrieve results strcuture from the workspace and check if it exists
disableFields = true;
resultStruct = getappdata(0,'optimTool_results_121677');
% If results structure does not exist then disable the 'Export results...' choice
if ~isempty(resultStruct)
    disableFields = false; 
end

% Call the function to layout the dialog box
[checkBoxes, editFields] = layoutDialog(hDialog, okButton, cancelButton, ...
                                        checkboxLabels, variableNames,disableFields);
% Set callback function for OK button
set(okButton, 'Callback', {@OKCallback, hDialog, checkBoxes, editFields, ...
                           optStruct, probStruct,resultStruct});
% Set callback function for keyboard responses
set(hDialog, 'KeyPressFcn', {@KeyPressCallback, hDialog, checkBoxes, editFields, ...
                           optStruct, probStruct,resultStruct});
%set the okButton to be the default
if (usejava('awt') == 1)
    fh = handle(hDialog);
    % Call the setDefaultButton method on the figure handle
    fh.setDefaultButton(okButton);
end

% Show dialog now!
set(hDialog, 'Visible', 'on');

% Set focus to the OK button
uicontrol(okButton);

%----------------------------------------------------------------------------
function modifiedNames = createVarNames(defVariableNames)
    % Preallocating for speed
    modifiedNames = cell(1, length(defVariableNames));
    for i = 1:length(defVariableNames)
        modifiedNames{i} = computename(defVariableNames{i});
    end

%----------------------------------------------------------------------------
function name = computename(nameprefix)

if (evalin('base',['exist(''', nameprefix,''', ''var'');']) == 0)
    name = nameprefix;
    return
end

% get all names that start with prefix in workspace
workvars = evalin('base', ['char(who(''',nameprefix,'*''))']);
% trim off prefix name
workvars = workvars(:,length(nameprefix)+1:end); 

if ~isempty(workvars)
    % remove all names with suffixes that are "non-numeric"
    lessthanzero = workvars < '0';
    morethannine = workvars > '9';
    notblank = (workvars ~= ' ');
    notnumrows = any((notblank & (lessthanzero | morethannine)),2);
    workvars(notnumrows,:) = [];
end

% find the "next one"
if isempty(workvars)
    name = [nameprefix, '1'];
else
    nextone = max(str2num(workvars)) + 1;
    if isempty(nextone)
        name = [nameprefix, '1'];
    else
        name = [nameprefix, num2str(nextone)];
    end
end

%----------------------------------------------------------------------------
function OKCallback(obj, eventdata, dialog, cb, e, optStruct, probStruct,resultStruct)

    CB_PROBLEM = 1;
    CB_RESTART = 2;
    CB_OPTION = 3;
    CB_RESULTS = 4;
    
    varnames = [];
    
     % we care only about items that are checked
     for i = 1:length(e)
         if get(cb{i}, 'Value') == 1 && i~=CB_RESTART
            varnames{end + 1} = get(e{i}, 'String');
         end
     end
    
     if isempty(varnames)
         errordlg(getString(message('optim:optimtool:ErrDlgSelectItemToExport')), ...
                  getString(message('optim:optimtool:ExportOptsDlgTitle')));
         return;
     end
    
    %check for invalid and empty variable names
    badnames = [];
    numbadentries = 0;
    emptystrmsg = '';
    badnamemsg = '';
    for i = 1:length(varnames)
        if strcmp('', varnames{i})
            numbadentries = numbadentries + 1;
            emptystrmsg = sprintf( ...
                getString(message('optim:optimtool:EmptyStrInvalidVarName')));
        elseif ~isvarname(varnames{i})
            badnames{end + 1} = varnames{i};
            numbadentries = numbadentries + 1;
        end
    end
    badnames = unique(badnames);
   
    if ~isempty(badnames)
        if (length(badnames) == 1)
            badnamemsg = getString(message('optim:optimtool:InvalidVarName',badnames{1}));
        elseif (length(badnames) == 2)
            badnamemsg = getString(message('optim:optimtool:TwoInvalidVarNames', ...
                      badnames{1},badnames{2}));
        else 
            badnamemsg = getString(message('optim:optimtool:MultInvalidVarNames', ...
                      sprintf('"%s", ',badnames{1:end-2}),badnames{end-1},badnames{end}));
        end
    end
    
    if numbadentries > 0 
        dialogname = getString(message('optim:optimtool:ErrDlgTitleInvalidVarNames'));
        if numbadentries == 1
            dialogname = getString(message('optim:optimtool:ErrDlgTitleInvalidVarName'));
        end
        errordlg([emptystrmsg badnamemsg], dialogname);    
        return; 
    end
    
    %check for names already in the workspace
    dupnames = [];
    for i = 1:length(varnames)
        if evalin('base',sprintf('exist(''%s'', ''var'');',varnames{i}))
            dupnames{end + 1} = varnames{i};
        end
    end
    dupnames = unique(dupnames);
 
    if ~isempty(dupnames) 
        dialogname = getString(message('optim:optimtool:ErrDlgTitleDupVarNames'));
        if (length(dupnames) == 1)
            queststr = getString(message('optim:optimtool:DupVarName',dupnames{1}));
            dialogname = getString(message('optim:optimtool:ErrDlgTitleDupVarName'));
        elseif (length(dupnames) == 2)
            queststr = getString(message('optim:optimtool:TwoDupVarNames', ...
                dupnames{1},dupnames{2}));
        else
            queststr = getString(message('optim:optimtool:MultDupVarNames', ...
                sprintf('"%s", ',dupnames{1:end-2}),dupnames{end-1},dupnames{end}));
        end
        buttonName = questdlg(queststr, dialogname,...
            getString(message('optim:optimtool:BtnYes')), ... % Yes
            getString(message('optim:optimtool:BtnNo')), ... % No
            getString(message('optim:optimtool:BtnYes')));  % default: Yes
        if ~strcmp(buttonName, getString(message('optim:optimtool:BtnYes'))) 
            return;
        end 
    end

    includeRestart = get(cb{CB_RESTART}, 'Value') == 1;
    %Check for variable names repeated in the dialog edit fields
    [uniqueArray ignore uniqueIndex] = unique(varnames);
    if length(varnames) == length(uniqueArray)
        if get(cb{CB_PROBLEM}, 'Value') == 1  % Export problem structure
            % Input argument 'probStruct' is modified so that it contains only fields relevant
            % to the solver 'probStruct.solver'
            probStruct = createProblemStruct(probStruct.solver,[],probStruct);
            probStruct.options = createOptionsStruct(probStruct.solver,optStruct);
            if includeRestart
                probStruct = addRestartInfo(probStruct.solver,probStruct,resultStruct);
            end
            assignin('base', get(e{CB_PROBLEM}, 'String'), probStruct);
        end
        if get(cb{CB_OPTION}, 'Value') == 1   % Export options structure
            optStruct = createOptionsStruct(probStruct.solver,optStruct);
            if includeRestart
                optStruct = addRestartInfo(probStruct.solver,optStruct,resultStruct);
            end
            assignin('base', get(e{CB_OPTION}, 'String'), optStruct);
        end
        if get(cb{CB_RESULTS}, 'Value') == 1  % export result structure
            assignin('base', get(e{CB_RESULTS}, 'String'), resultStruct); 
        end
        if length(varnames) == 1
            msg = getString(message('optim:optimtool:CreatedVarInWS',varnames{1}));
        elseif length(varnames) == 2 
            msg = getString(message('optim:optimtool:CreatedVarsInWS',varnames{1},varnames{2}));
        elseif length(varnames) == 3
            msg = getString(message('optim:optimtool:CreatedMultVarsInWS',varnames{1},varnames{2},varnames{3}));
        else  %shouldn't get here
            msg='';
        end
        disp(msg);
        delete(dialog);
    else
        errordlg(getString(message('optim:optimtool:NotUniqueNames')),...
            getString(message('optim:optimtool:ErrDlgTitleNotUniqueNames')));
    end
%----------------------------------------------------------------------------
function myStruct = addRestartInfo(solver,myStruct,resultStruct)
solversWithX0 = {'fmincon','fminunc','lsqnonlin', ...
    'lsqcurvefit','linprog','quadprog', ...
    'bintprog','fgoalattain','fminimax', ...
    'fseminf','fminsearch','fzero', ...
    'fsolve','lsqlin','lsqnonneg', ...
    'patternsearch','simuannealbnd'};
solversWithPop = {'ga','gamultiobj'};
% fminbnd does not have any problem data that can help restart from the
% final solution
if any(strcmpi(solver,solversWithX0)) && isfield(myStruct,'solver')
    myStruct.x0 = resultStruct.x;
end
if any(strcmpi(solver,solversWithPop))
    if isfield(myStruct,'solver')
        myStruct.options.InitialPopulation = resultStruct.population;
        myStruct.options.InitialScores = resultStruct.score;
    else
        myStruct.InitialPopulation = resultStruct.population;
        myStruct.InitialScores = resultStruct.score;
    end
end

%----------------------------------------------------------------------------
function CancelCallback(obj, eventdata, dialog)
    delete(dialog);
   
%----------------------------------------------------------------------------
function KeyPressCallback(obj, eventdata, dialog, cb, e, optStruct, probStruct,resultStruct)
% This function is a callback for response from keyboard instead of mouse (a wrapper around 
% 'OKCallback' function)
    asciiVal = get(dialog, 'CurrentCharacter');
    if ~isempty(asciiVal)
        % space bar or return is the "same" as OK
        if (asciiVal==32 || asciiVal==13)   
            OKCallback(obj, eventdata, dialog, cb, e, optStruct, probStruct,resultStruct);
        elseif (asciiVal == 27) % Escape has the same effect as Cancel
            delete(dialog);
        end
    end
   
%----------------------------------------------------------------------------
function [cb, e] = layoutDialog(hDlg, okBut, cancelBut, checkboxLabels, ...
                                variableNames,disableFields)
% Dialog position and other properties are set in this function

    EXTENT_WIDTH_INDEX = 3;  % width is the third argument of extent
    
    POS_X_INDEX      = 1;
    POS_Y_INDEX      = 2;
    POS_WIDTH_INDEX  = 3;
    POS_HEIGHT_INDEX = 4;
    
    CONTROL_SPACING  = 10;
    EDIT_WIDTH       = 90;
    CHECK_BOX_WIDTH  = 20;
    DEFAULT_INDENT   = 20;
    
    numItem = 4;
    CB_PROBLEM = 1;
    CB_RESTART = 2;
    CB_OPTION = 3;
    CB_RESULTS = 4;
    
    % the following mimics what questdlg does in terms of setting the
    % height of the buttons, which are a little larger than the default
    % uicontol button height
    
    % questdlg also sets the y position of the buttons to 10.
    
    BtnYOffset = 10;
    
    BtnHeight = 22;
    BtnMargin = 1.4;
    BtnExtent = get(okBut, 'extent');
    BtnHeight = max(BtnHeight,BtnExtent(4)*BtnMargin);
    
    okPos = get(okBut, 'Position');
    okPos(POS_HEIGHT_INDEX) = BtnHeight;
    okPos(POS_Y_INDEX) = BtnYOffset;
    set(okBut, 'Position', okPos);
    cancelPos = get(cancelBut, 'Position');
    cancelPos(POS_HEIGHT_INDEX) = BtnHeight;
    cancelPos(POS_Y_INDEX) = BtnYOffset;
    set(cancelBut, 'Position', okPos);
    
    longestCBExtent = 0;
    ypos = okPos(POS_HEIGHT_INDEX) + okPos(POS_Y_INDEX)+ 2*CONTROL_SPACING;
    cb = cell(numItem, 1);
    e = cell(numItem, 1);
    cbTag = {'ProbsAndOptsChkBox','InfoToResumeChkBox','OptsChkBox',...
             'ResultsChkBox'};
    ebTag = {'ProbsAndOptsEdtBox','UnusedEdtBox','OptsEdtBox',...
             'ResultsEdtBox'};
    for i = numItem:-1:1
        cb{i} = uicontrol(hDlg, 'Style', 'checkbox', 'String', ...
                          checkboxLabels{i}, 'Unit', 'pixels', ...
                          'Tag',cbTag{i});
        check_pos = get(cb{i}, 'Position');
        check_pos(POS_Y_INDEX) = ypos;
        extent = get(cb{i}, 'Extent');
        width = extent(EXTENT_WIDTH_INDEX);
        check_pos(POS_WIDTH_INDEX) = width + CHECK_BOX_WIDTH;  
        set(cb{i}, 'Position', check_pos);
        
        if (i==CB_RESTART) % indent 2nd line a little and don't add a edit field;
            check_pos(POS_X_INDEX) = check_pos(POS_X_INDEX) + CHECK_BOX_WIDTH;
            set(cb{i}, 'Position', check_pos);
        else
            e{i} = uicontrol(hDlg, 'Style', 'edit', 'String', variableNames{i}, ...
                'BackgroundColor', 'white', 'Unit', 'pixels', ...
                'HorizontalAlignment', 'left', 'Tag', ebTag{i});

            edit_pos = get(e{i}, 'Position');
            edit_pos(POS_Y_INDEX) = ypos;
            edit_pos(POS_WIDTH_INDEX) = EDIT_WIDTH;
            % cursor doesn't seem to appear in default edit height
            edit_pos(POS_HEIGHT_INDEX) = edit_pos(POS_HEIGHT_INDEX) + 1;
            set(e{i}, 'Position', edit_pos);
        end
        ypos = ypos + CONTROL_SPACING + edit_pos(POS_HEIGHT_INDEX);
        if width > longestCBExtent
            longestCBExtent = width;
        end
        
        if disableFields
            set(cb{CB_RESTART}, 'Enable', 'off');
            set(e{CB_RESTART}, 'Enable', 'off');
            set(e{CB_RESTART}, 'Backgroundcolor', [0.831373 0.815686 0.784314]);            
            set(cb{CB_RESULTS}, 'Enable', 'off');
            set(e{CB_RESULTS}, 'Enable', 'off');
            set(e{CB_RESULTS}, 'Backgroundcolor', [0.831373 0.815686 0.784314]);
            if strcmp(get(cb{CB_PROBLEM}, 'Enable'), 'off')  % only options is enabled - check it
                set(cb{CB_OPTION}, 'Value', 1);
            end
        end
    end

    % Position edit boxes
    edit_x_pos = check_pos(POS_X_INDEX) + longestCBExtent + CONTROL_SPACING ...
                           + CHECK_BOX_WIDTH;
    for i = 1:numItem
        edit_pos = get(e{i}, 'Position');
        edit_pos(POS_X_INDEX) = edit_x_pos;
        set(e{i}, 'Position', edit_pos);
    end
    h_pos = get(hDlg, 'Position');
    
    h_pos(POS_WIDTH_INDEX) = max(edit_x_pos + edit_pos(POS_WIDTH_INDEX) + ...
                                 CHECK_BOX_WIDTH, okPos(POS_WIDTH_INDEX) + ...
                                 cancelPos(POS_WIDTH_INDEX) + ...
                                 CONTROL_SPACING + (2 * DEFAULT_INDENT));
    h_pos(POS_HEIGHT_INDEX) = ypos;
    
    % Center it
    oldu = get(0,'Units');
    set(0,'Units','pixels');
    screenSize = get(0,'ScreenSize');
    set(0,'Units',oldu);
    
    h_pos(POS_X_INDEX) = (screenSize(3) - h_pos(3))/2;
    h_pos(POS_Y_INDEX) = (screenSize(4) - h_pos(4))/2;
    
    set(hDlg, 'Position', h_pos);
     
    % Make sure it is on-screen
    outerPos = get(hDlg,'OuterPosition');
    if outerPos(1)+outerPos(3) > screenSize(3)
        outerPos(1) = screenSize(3) - outerPos(3);
    end
    if outerPos(2)+outerPos(4) > screenSize(4)
        outerPos(2) = screenSize(4) - outerPos(4);
    end
    set(hDlg, 'OuterPosition', outerPos);
    
    x_ok = (h_pos(POS_WIDTH_INDEX))/2 -  (okPos(POS_WIDTH_INDEX) + ... 
            CONTROL_SPACING + cancelPos(POS_WIDTH_INDEX))/2;
    okPos(POS_X_INDEX) = x_ok;
    set(okBut, 'Position', okPos);
    cancelPos(POS_X_INDEX) = okPos(POS_X_INDEX) + okPos(POS_WIDTH_INDEX) + ...
                                   CONTROL_SPACING;
    set(cancelBut, 'Position', cancelPos);

    % Reorder the children so that tabbing makes sense
    children = get(hDlg, 'children');
    children = flipud(children);
    set(hDlg, 'children', children);

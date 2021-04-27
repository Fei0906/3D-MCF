function stop = optimtooloutput(~,optimval,state,varargin)
%

%OPTIMTOOLOUTPUT OutputFcn used to interact between the GUI and solvers.
%  It is used to capture and react to changes in the state of the GUI 
%  (STOP, PAUSE, RESUME, KILL). It also manages warnings by displaying them
%  in Status and Results panel. It also sets the iteration number.
%
%   Private to OPTIMTOOL.

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/06/30 16:36:57 $

%Initialize
stop = false;
drawnow;
% Three buttons in the GUI
STOP = 0;
RUN_RESUME = 1;
PAUSE = 2;

optimtoolGui = javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'); % Get a handle to the GUI
if isempty(optimtoolGui)
    stop = true;
    return;
end

switch state
    case 'init'
        [msg,id] = lastwarn;
        if ~isempty(msg)
            % Remove html tags from the warning message
            msg = regexprep(msg,'</?(\w+).*?>','');
            % Make sure that exactly one newline character is at the end of 'msg'
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            javaMethodEDT('appendResults',optimtoolGui, ...   % Append to the 'Status and Results' panel
                [getString(message('optim:optimtool:OptimtoolStatusWarning')),' ',msg]); 
            msg(regexp(msg,'\n')) = ' '; % Remove all newline characters before showing 'msg' in the dialog
            warndlg(msg,'Optimization Tool'); % Show a warning dialog box
        end
        % To avoid displaying warnings with same ID multiple times we store the warning IDs in appdata.
        % We show it only once in the 'Status and Results' panel.
        setappdata(0,'last_warning_id_for_optimtool',id);
        return; % Nothing to do now
    case {'iter', 'interrupt'}
        [msg,id] = lastwarn;
        % In addition to non-empty message the warning ID must be new (not in the appdata)
        if ~isempty(msg) && ~strcmp(id,getappdata(0,'last_warning_id_for_optimtool'))
            % Remove html tags from the warning message
            msg = regexprep(msg,'</?(\w+).*?>','');
            % Add a newline to the message if it does not have already
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            javaMethodEDT('appendResults',optimtoolGui, ...
                [getString(message('optim:optimtool:OptimtoolStatusWarning')),' ',msg]);
            lastwarn('');
        end
        javaMethodEDT('setIteration',optimtoolGui,value2RHS(optimval.iteration)); % Update iteration number in the GUI
        % Action based on run mode of GUI
        RunMode = javaMethodEDT('getRunMode',optimtoolGui);
        switch RunMode
            case RUN_RESUME

            case STOP
                stop = true;
                return;
            case PAUSE
                fprintf(getString(message('optim:optimtool:OptimtoolPaused')));
                % If in pause state keeping looping here.
                while true
                    drawnow
                    if isempty(javaMethodEDT('getOptimGUI','com.mathworks.toolbox.optim.OptimGUI'))
                        stop = true; % The GUI was closed by user
                        return;
                    end
                    mode = javaMethodEDT('getRunMode',optimtoolGui);
                    if mode == STOP
                        stop = true;
                        return;
                    elseif mode == RUN_RESUME
                        break;
                    end
                end % End while
            otherwise % Safeguard
                return;
        end
    case 'done'
        javaMethodEDT('setIteration',optimtoolGui,value2RHS(optimval.iteration)); % Update iteration number in the GUI
        if isappdata(0,'last_warning_id_for_optimtool')
            rmappdata(0,'last_warning_id_for_optimtool');
        end
        return;
    otherwise % Safeguard
        return;
end

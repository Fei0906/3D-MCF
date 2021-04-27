function hashModel = setoptimrandstates(problemStruct,hashModel,randchoice)
%SETOPTIMRANDSTATES set states of random number generator for solvers
%
%   Private to OPTIMTOOL

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/05/10 17:32:12 $

% Get structure that stores the states
if isappdata(0,'optim_rand_states_struct')
    optimSolversRandStates = getappdata(0,'optim_rand_states_struct');
else % Create a structure that can store all the random states
    optimSolversRandStates.ga.rngstate = [];
    optimSolversRandStates.ga.garandchoice = false;

    optimSolversRandStates.gamultiobj.rngstate = [];
    optimSolversRandStates.gamultiobj.gamultiobjrandchoice = false;

    optimSolversRandStates.patternsearch.rngstate = [];
    optimSolversRandStates.patternsearch.patternsearchrandchoice = false;

    optimSolversRandStates.simulannealbnd.rngstate = [];
    optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = false;

    % Save in the appdata
    setappdata(0,'optim_rand_states_struct',optimSolversRandStates);
end

% Return if the hashModel does not have 'rngstate' key
% Any modification in the hashModel is okay when nargout == 0; In this
% case, simply set the states in the appdata for the correct solver
if nargout > 0 && ~hashModel.containsKey('rngstate')
    return;
else
    % Remove rngstate from the hashModel. The GUI only
    % needs to know if rngstate is present or not (a boolean flag
    % 'randchoice') 
    hashModel.remove('rngstate');
    if randchoice
        useCurrentRandState = true;
    else
        useCurrentRandState = false;
    end
end

% Get the random states from the problem structure and save in the appdata
rngstate = problemStruct.rngstate;

if ~isempty(rngstate)
    switch problemStruct.solver
        case 'ga'
            optimSolversRandStates.ga.rngstate = rngstate;
        case 'gamultiobj'
            optimSolversRandStates.gamultiobj.rngstate = rngstate;
        case 'patternsearch'
            optimSolversRandStates.patternsearch.rngstate = rngstate;
        case 'simulannealbnd'
            optimSolversRandStates.simulannealbnd.rngstate = rngstate;
    end
end

switch problemStruct.solver
    case 'ga'
        if hashModel.containsKey('garandchoice') 
            if strcmpi('true',hashModel.get('garandchoice'))
                optimSolversRandStates.ga.garandchoice = true;
            elseif strcmpi('false',hashModel.get('garandchoice'))
                optimSolversRandStates.ga.garandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.ga.garandchoice = true;
            hashModel.put('garandchoice','true');
        end
    case 'gamultiobj'
        if hashModel.containsKey('gamultiobjrandchoice')
            if strcmpi('true',hashModel.get('gamultiobjrandchoice'))
                optimSolversRandStates.gamultiobj.gamultiobjrandchoice = true;
            elseif strcmpi('false',hashModel.get('gamultiobjrandchoice'))
                optimSolversRandStates.gamultiobj.gamultiobjrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.gamultiobj.gamultiobjrandchoice = true;
            hashModel.put('gamultiobjrandchoice','true');
        end
    case 'patternsearch'
        if hashModel.containsKey('patternsearchrandchoice')
            if strcmpi('true',hashModel.get('patternsearchrandchoice'))
                optimSolversRandStates.patternsearch.patternsearchrandchoice = true;
            elseif strcmpi('false',hashModel.get('patternsearchrandchoice'))
                optimSolversRandStates.patternsearch.patternsearchrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.patternsearch.patternsearchrandchoice = true;
            hashModel.put('patternsearchrandchoice','true');
        end
    case 'simulannealbnd'
        if hashModel.containsKey('simulannealbndrandchoice') 
            if strcmpi('true',hashModel.get('simulannealbndrandchoice'))
                optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = true;
            elseif strcmpi('false',hashModel.get('simulannealbndrandchoice'))
                optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = false;
            end
        elseif useCurrentRandState
            optimSolversRandStates.simulannealbnd.simulannealbndrandchoice = true;
            hashModel.put('simulannealbndrandchoice','true');
        end
end
% Save in the appdata
setappdata(0,'optim_rand_states_struct',optimSolversRandStates);

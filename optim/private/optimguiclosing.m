function optimguiclosing()
%OPTIMGUICLOSING GUI helper function to clean up 'result' appdata from 
%   the MATLAB workspace when the GUI is closed.

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:50:39 $

% Remove appdata structures used by optimtool
if isappdata(0,'optimTool_results_121677')
    rmappdata(0,'optimTool_results_121677');
end

if isappdata(0,'optimTool_Problem_Data')
    rmappdata(0,'optimTool_Problem_Data');
end

if isappdata(0,'optimTool_Options_Data')
    rmappdata(0,'optimTool_Options_Data');
end

if isappdata(0,'optimTool_Problem_HashTable')
    rmappdata(0,'optimTool_Problem_HashTable');
end

if isappdata(0,'optimTool_Options_HashTable')
    rmappdata(0,'optimTool_Options_HashTable');
end

if isappdata(0,'optim_rand_states_struct')
    rmappdata(0,'optim_rand_states_struct');
end

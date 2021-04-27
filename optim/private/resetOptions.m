function resetOptions()
% RESETOPTIONS Reset OPTIONS
%   Private to OPTIMTOOL

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/12/10 21:50:45 $

% Create the default options structure
options = createOptionsStruct('all',struct('Display', 'off')); % This is the default for the GUI
setappdata(0,'optimTool_Options_Data',options);
resetOptimtoolHashTable('optimTool_Options_HashTable');



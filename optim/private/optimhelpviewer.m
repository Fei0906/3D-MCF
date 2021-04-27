function optimhelpviewer()
% OPTIMHELPVIEWER is a helper file for Optimtool. 

%   Copyright 2006-2011 The MathWorks, Inc. 
%   $Revision: 1.1.6.4 $

mapfilename = fullfile(docroot,'toolbox','optim','helptargets.map');
try
    helpview(mapfilename, 'optimtool');
catch
    msg = getString(message('optim:optimtool:UnableToOpenHelp'));
    errordlg(msg);
end

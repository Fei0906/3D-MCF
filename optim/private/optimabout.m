function optimabout()
%OPTIMABOUT helper that displays the About Box  

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/06/30 16:37:05 $

a = ver('optim');
str = sprintf(['Optimization Toolbox %s\n',...
               'Copyright 1990-%s The MathWorks, Inc.'], ...
               a.Version,a.Date(end-3:end));
aboutTitle = getString(message('optim:optimtool:TitleAboutOptimTlbx'));
msgbox(str,aboutTitle,'modal');

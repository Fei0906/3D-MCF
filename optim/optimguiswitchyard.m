function varargout = optimguiswitchyard(action,varargin)
%

%OPTIMGUISWITCHYARD switchyard for optimtool.
% Helper function for the optimtool

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ $Date: 2011/05/09 01:06:00 $

% Calls from Java prefer the if/else version.
% [varargout{1:max(nargout,1)}]=feval(action,varargin{:});
if nargout==0
	feval(action,varargin{:});
else    
	[varargout{1:nargout}]=feval(action,varargin{:});
end

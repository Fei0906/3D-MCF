function relDelta = relDeltaXFminunc(x,deltaX)
%relDeltaXFminunc Helper function for fminunc.
% 
% Calculates the relative displacement, as used in fminunc medium-scale
% stopping test.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/08/29 20:34:17 $

relDelta = norm(deltaX ./ (1 + abs(x)),inf);

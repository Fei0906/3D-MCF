function displayProblemInfo(sizes)
% diagnose is a helper function that displays diagnostics information

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2008/09/13 06:57:34 $

fprintf('\n')
fprintf(' Variables:                        \t%10i\n',sizes.nVar)
fprintf(' Finite lower bounds:              \t%10i\n',sizes.nFiniteLb)
fprintf(' Finite upper bounds:              \t%10i\n',sizes.nFiniteUb)
fprintf(' Fixed variables:                  \t%10i\n',sizes.nFixedVar)
fprintf('\n')
fprintf(' Linear equality constraints:      \t%10i\n',sizes.mLinEq)
fprintf(' Linear inequality constraints:    \t%10i\n',sizes.mLinIneq)
fprintf('\n')
fprintf(' Nonlinear equality constraints:   \t%10i\n',sizes.mNonlinEq)
fprintf(' Nonlinear inequality constraints: \t%10i\n',sizes.mNonlinIneq)
fprintf('\n')

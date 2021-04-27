function [cModified, AModified, bModified, lbModified, ubModified, xModified, exitflag] = ...
         changeBoundDualFeasible(chglbidx, chglbval, chgubidx, chgubval, c, A, b,... 
              Aeq, beq, lb, ub, x, basicVarIdx, nonbasicVarIdx, delrows, tol)
%changeBoundDualFeasible Process bounds and construct dual feasible point. 
% This function changes the lower and upper bounds of the original problem 
% and construct a dual feasible point to the modified problem with changed bounds.              
% min  c'*x 
% s.t.  A x <= b
%       Aeq x = beq
%       lb <= x <= ub
%
% EXITFLAG describes the exit conditions:
%   If EXITFLAG is: 
%      = 1 the changed bounds are acceptable and the reconstructed point is
%      dual feasible to the modified system with only equality constraints.
%      = -1 the changed bounds were not compatible and the modified problem
%      became infeasible.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2010/10/08 17:13:33 $

% Initialization.
exitflag = 1;
mineq = size(A, 1);
meq = size(Aeq, 1);

% Transform into the standard form with only equality constraints and bounds. 
% The standard form is
% min  cModified'* y 
% s.t.  AModified * y = bModified
%       lbModified <= y <= ubModified.
% Add slack variables, add their bounds and objective coefficients accordingly.
AModified = [A  eye(mineq); Aeq  zeros(meq,mineq)];
bModified = [b; beq];
cModified = [c; zeros(mineq,1)];
lbModified = [lb; zeros(mineq,1)];
ubModified = [ub; inf(mineq, 1)];

% Construct a point xModified = [x; b-A*x] according to the new system.
if ~isempty(A)
    xModified = [x; b-A*x];
else
    xModified = x;
end
% Assumption: the input basicVarIdx and nonbasicVarIdx include the indices 
% of slack variables.

% Check that chglbidx <= nvars and chgubidx <= nvars
% Check the changed bound satisfies lbModified <= ubModified.
nvars  = length(c);
if any(chglbidx > nvars) || any(chgubidx >nvars)
    error(message('optim:changeBoundDualFeasible:IdxExceedNumofVars'));
end
% Change of bounds.
lbModified(chglbidx) = chglbval;
ubModified(chgubidx) = chgubval;

if any(lbModified > ubModified)
    error(message('optim:changeBoundDualFeasible:IncompatibleBounds'));
end

if any(delrows)
    leftrows = ~delrows; 
    AModified = AModified(leftrows, :);
    bModified = bModified(leftrows);
end

% Construct a new dual feasible basis as the starting basis for dual simplex phase two.
% Reset the nonbasic variables (only those with changed bounds).
tol2 = tol^2;
xAtLower = abs(x(chglbidx) - lb(chglbidx)) <= tol2;
xAtUpper = abs(x(chgubidx) - ub(chgubidx)) <= tol2;
% Only needed for nonbasic variables, but ok for basic since basic variables 
% will be overwritten.
xModified(chglbidx(xAtLower)) = chglbval(xAtLower);
xModified(chgubidx(xAtUpper)) = chgubval(xAtUpper);

% If xModified(i)=lb(i)=ub(i)), reset xModified(i) at the new lower or upper bounds 
% depending on the sign of dualf(i) so that xModified will be a starting dual feasibility point.
if any( (ub-lb) <= tol2)
    y0 = cModified(basicVarIdx)' / AModified(:, basicVarIdx);
    dualf0 = cModified(nonbasicVarIdx) - (y0*AModified(:, nonbasicVarIdx))'; % dual feasibility for the original system.
    dualf0p = dualf0 > tol;
    dualf0n = dualf0 < -tol;
    xModified(nonbasicVarIdx(dualf0p)) = lbModified(nonbasicVarIdx(dualf0p));
    xModified(nonbasicVarIdx(dualf0n)) = ubModified(nonbasicVarIdx(dualf0n));
end

B = AModified(:, basicVarIdx);
N = AModified(:, nonbasicVarIdx);
if ~isempty(B)
    xModified(basicVarIdx) = B \ (bModified - N* xModified(nonbasicVarIdx));
end

if ~isempty(bModified) && max(abs(AModified*xModified - bModified)) > tol
    error(message('optim:changeBoundDualFeasible:XModifiedViolatesConstr'));
end
%--------------------------------------------------------------------------
%--------------END of changeBoundDualFeasible------------------------------
%--------------------------------------------------------------------------
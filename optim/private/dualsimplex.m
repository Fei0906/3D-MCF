function [xsol, fval, exitflag, output, basicVarIdx, nonbasicVarIdx] = ...
    dualsimplex(c, Aeq, beq, lb, ub, x0, basicVarIdx0, nonbasicVarIdx0, maxiter, tol, verbosity)
% DUALSIMPLEX  Phase two of dual simplex method for Linear Programming
%               min   c' * x 
%        subject to  Aeq * x  = beq; 
%                    lb <= x <= ub.
%
%   It starts with a dual feasible point x0 with its corresponding basic
%   variable index basicVarIdx0 and nonbasic variable index nonbasicVarIdx0
%   under the assumption that the above standard form of linear programming is used.
%
%   EXITFLAG describes the exit conditions:
%   If EXITFLAG is:
%      = 1 then DUALSIMPLEX converged with a solution XSOL.
%      = 0 then DUALSIMPLEX did not converge within the maximum number of
%      iterations.
%      = -1 then the problem was infeasible.
%      = -2 then the problem was unbounded.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision.1 $  $Date: 2010/10/08 17:13:39 $

if isempty(x0) 
    error(message('optim:dualsimplex:EmptyStartPoint')); 
end

% Initialize the output arguments.
output.iterations = 0;
output.cgiterations = [];
output.algorithm = 'medium-scale: dual simplex';
basicVarIdx = basicVarIdx0(:);
nonbasicVarIdx = nonbasicVarIdx0(:);

% Initialization of constants for exitflag
Converged   = 1;
ExcdMaxiter = 0;
Infeasible  = -1;
Unbounded   = -2;
Unset       = inf;

% Initialization: a dual feasible basis
x_B = x0(basicVarIdx);
x_N = x0(nonbasicVarIdx);
lb_B = lb(basicVarIdx);
ub_B = ub(basicVarIdx);
lb_N = lb(nonbasicVarIdx);
ub_N = ub(nonbasicVarIdx);
c_B = c(basicVarIdx);
c_N = c(nonbasicVarIdx);
B = Aeq(:, basicVarIdx);
N = Aeq(:, nonbasicVarIdx);

% Initialization of dualf, the coefficients for nonbasic variables.
y = c_B' / B;
dualf = c_N - (y*N)'; 
% Check that x0 is a dual feasible basis solution, satisfying
% if dualf(j) > 0, then x_N(j) == lb_N(j)
% if dualf(j) < 0, then x_N(j) == ub_N(j)
checkconstr = all(Aeq *x0 - beq <= tol);
checkdualfeas =  all( (abs(dualf) < tol) | ((dualf >= -tol) & ((x_N - lb_N) <= tol)) ...
                               | ((dualf <= tol) & ((ub_N - x_N) <= tol)) );
if ~(checkconstr && checkdualfeas)
     error(message('optim:dualsimplex:DualInfeasibleStartPoint'));
end

% Initialize output variables.
niters = 0;
fval = c'*x0;
xsol = [x_B; x_N];
exitflag = Unset;

% Display column headers.
if verbosity >= 2
    disp( sprintf('\nMinimize using dual simplex.') );
    disp( sprintf('                                             Maximal Primal') );
    disp( sprintf('      Iter            Objective              Infeasibility') );
    rbounds = max(lb_B - x_B, 0) + max(x_B - ub_B, 0); 
    disp( sprintf('%8d         %12.6g                %12.6g', niters, fval, norm(rbounds, inf)) );
end

%------------------------ while loop ---------------------
while (niters < maxiter) 
    
    % Check the optimality conditions: primal feasibility.
    pfeas = (x_B <= ub_B + tol) & (x_B >= lb_B - tol);
    if all( pfeas )
        xsol = [x_B; x_N]; % xsol should have already been computed.
        fval = [c_B; c_N]'* xsol;
        exitflag = Converged;
        if verbosity > 0
            disp('Optimization terminated.');
        end
        break;
    end
    
    % Otherwise, choose a leaving variable.
    indexLeave = find((~pfeas), 1, 'first');  
    
    % Solve the linear system to choose the entering variable.
    rhs = zeros(length(basicVarIdx),1); 
    rhs(indexLeave) = 1;
    v = rhs'/B; % Note v is a row vector.
    w = (v * N)';
    
    % Collect all the candidates for entering variables.
    if x_B(indexLeave) < lb_B(indexLeave)
        entCands = ((w < -tol) & (x_N < ub_N -tol)) | ((w > tol) & (x_N > lb_N + tol)); 
    elseif x_B(indexLeave) > ub_B(indexLeave)
        entCands = ((w > tol) & (x_N < ub_N - tol)) | ((w < -tol) & (x_N > lb_N + tol));
    else
        error(message('optim:dualsimplex:WrongChoiceOfLeavingVar'));
    end
    
    % Choose an entering variable.
    if ~any(entCands)
        % If no entering candidates exist, the problem is infeasible.
        exitflag = Infeasible ;
        if verbosity > 0
            disp('Exiting: The constraints are overly stringent; no feasible starting point found.');
        end
        break;
    else
        % Otherwise, choose the index that minimizes |dualf_j/ w_j|,
        % where j is the index of an entering variable candidate.
        [tmp, tindx] = min( abs( dualf(entCands)./w(entCands) ) );
        eindx = find(entCands);
        indexEnter = eindx(tindx); % the index of nonbasicVarIdx        
    end
    
    % Solve the linear system for updating x.
    N_indexEnter = N(:, indexEnter);
    d = B\N_indexEnter;
    
    % Compute the value of delta in two cases:
    if x_B(indexLeave) < lb_B(indexLeave)
        delta = (x_B(indexLeave) - lb_B(indexLeave))/w(indexEnter);
    elseif x_B(indexLeave) > ub_B(indexLeave) 
        delta = (x_B(indexLeave) - ub_B(indexLeave))/w(indexEnter);
    end
       
    % Update the dual feasible basis solution
    x_N(indexEnter) = x_N(indexEnter) + delta;
    x_B = x_B - delta * d;

    %----------------------------------------------------------------------
    % Update Basis, Nonbasis; basicVarIdx, nonbasicVarIdx; and c_B, c_N; 
    % ub_B, ub_N; lb_B, lb_N; 
    % x_B(indexLeave), x_N(indexEnter) after each iteration on basis
    % Basis(:,indexLeave) <--> Nonbasis(:, indexEnter);
    % basicVarIdx(indexLeave)   <--> nonbasicVarIdx(indexEnter);
    % c_B(indexLeave)     <--> c_N(indexEnter);
    % ub_B(indexLeave)    <--> ub_N(indexEnter);
    % lb_B(indexLeave)    <--> lb_N(indexEnter);
    % x_B(indexLeave)     <--> x_N(indexEnter);
    %---------------------------------------------------------------------
    N(:, indexEnter) = B(:, indexLeave); 
    B(:, indexLeave) = N_indexEnter;
    
    swaptmp = basicVarIdx(indexLeave);
    basicVarIdx(indexLeave) = nonbasicVarIdx(indexEnter);
    nonbasicVarIdx(indexEnter) = swaptmp;
    
    swaptmp = x_B(indexLeave);
    x_B(indexLeave) = x_N(indexEnter);
    x_N(indexEnter) = swaptmp;
    
    swaptmp = lb_B(indexLeave);
    lb_B(indexLeave) = lb_N(indexEnter);
    lb_N(indexEnter) = swaptmp;
    
    swaptmp = ub_B(indexLeave);
    ub_B(indexLeave) = ub_N(indexEnter);
    ub_N(indexEnter) = swaptmp;
    
    swaptmp = c_B(indexLeave);
    c_B(indexLeave) = c_N(indexEnter);
    c_N(indexEnter) = swaptmp;
     
    % Reset the value of dualf.
    t = - dualf(indexEnter)/w(indexEnter);
    dualf = dualf + t * w; 
    dualf(indexEnter) = t;
   
    % Compute and display the objective function value at each iteration.
    niters = niters + 1;
    xsol = [x_B; x_N];
    fval = [c_B; c_N]'*xsol;
    
    if verbosity >= 2
        rbounds = max(lb_B - x_B, 0) + max(x_B - ub_B, 0); 
        disp( sprintf('%8d         %12.6g                %12.6g', ...
                      niters, full(fval), norm(rbounds, inf)) );
    end
    
end % while (niters <= maxiter)

% Assignments to output arguments
[tmp, indx] = sort([ basicVarIdx; nonbasicVarIdx]);
xsol = xsol(indx);
output.iterations = niters;

if (niters == maxiter) && (exitflag == Unset)
    exitflag = ExcdMaxiter;
    if verbosity > 0
        disp('Exiting: Maximum number of iterations exceeded; increase options.MaxIter.');
    end
end

% Check the boundary conditions and eliminate roundoff error.
if (exitflag == Converged)
    closelb = (abs(xsol - lb) <= eps(max(abs(lb), 1))*100) & (xsol ~= lb);
    closeub = (abs(xsol - ub) <= eps(max(abs(ub), 1))*100) & (xsol ~= ub);
    xcloselb = xsol(closelb);
    xcloseub = xsol(closeub);
    if ~isempty([closelb; closeub])     
        xsol(closelb) = lb(closelb);
        xsol(closeub) = ub(closeub);
        
        checkAeq = abs(Aeq *xsol - beq) <= tol;
        if all(checkAeq)
           fval = c'*xsol;
        else %Return the computed solution.
           xsol(closelb) = xcloselb;
           xsol(closeub) = xcloseub;
        end
    end
end 
%--------------------------------------------------------------------------
%-------------END of DUALSIMPLEX ------------------------------------------
%--------------------------------------------------------------------------

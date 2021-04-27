function [x,fval,exitflag,output] = bintprog(f,A,b,Aeq,beq,x0,options)
%BINTPROG Binary integer programming.
%   BINTPROG solves the binary integer programming problem
%
%       min f'*X  subject to:  A*X <= b,
%                              Aeq*X = beq,
%                              where the elements of X are binary     
%                              integers, i.e., 0's or 1's.             
%
%   X = BINTPROG(f) solves the problem min f'*X, where the elements of X 
%   are binary integers. 
%
%   X = BINTPROG(f,A,b) solves the problem min f'*X subject to the linear 
%   inequalities A*X <= b, where the elements of X are binary integers.
%
%   X = BINTPROG(f,A,b,Aeq,beq) solves the problem min f'*X subject to the
%   linear equalities Aeq*X = beq, the linear inequalities A*X <= b, where 
%   the elements of X are binary integers.
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0) sets the starting point to X0. The 
%   starting point X0 must be binary integer and feasible, or it will 
%   be ignored.
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0,OPTIONS) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS, an
%   argument created with the OPTIMSET function.  See OPTIMSET for details.
%   Available options are BranchStrategy, Diagnostics, Display, 
%   NodeDisplayInterval, MaxIter, MaxNodes, MaxRLPIter, MaxTime, 
%   NodeSearchStrategy, TolFun, TolXInteger and TolRLPFun.
%
%   X = BINTPROG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the vector 'f' in PROBLEM.f, the linear inequality
%   constraints in PROBLEM.Aineq and PROBLEM.bineq, the linear equality
%   constraints in PROBLEM.Aeq and PROBLEM.beq, the start point in 
%   PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'bintprog' in PROBLEM.solver. Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. The structure PROBLEM 
%   must have all the fields. 
%
%   [X,FVAL] = BINTPROG(...) returns the value of the objective function at
%   X: FVAL = f'*X.
%
%   [X,FVAL,EXITFLAG] = BINTPROG(...) returns an EXITFLAG that describes 
%   the exit condition of BINTPROG. Possible values of EXITFLAG and the 
%   corresponding exit conditions are
%
%      1  BINTPROG converged to a solution X.
%      0  Maximum number of iterations exceeded.
%     -2  Problem is infeasible.
%     -4  MaxNodes reached without converging.              
%     -5  MaxTime reached without converging.               
%     -6  Number of iterations performed at a node to solve the 
%          LP-relaxation problem exceeded MaxRLPIter reached without 
%          converging.            
%
%   [X,FVAL,EXITFLAG,OUTPUT] = BINTPROG(...) returns a structure OUTPUT 
%   with the number of iterations in OUTPUT.iterations, the number of 
%   nodes explored in OUTPUT.nodes, the execution time (in seconds) in 
%   OUTPUT.time, the algorithm used in OUTPUT.algorithm, the branch 
%   strategy in OUTPUT.branchStrategy, the node search strategy in 
%   OUTPUT.nodeSrchStrategy, and the exit message in OUTPUT.message.
%
%   Example
%     f = [-9; -5; -6; -4]; 
%     A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
%     b = [9; 1; 0; 0];
%     X = bintprog(f,A,b) 
%
% See also LINPROG.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision.4 $  $Date: 2011/06/30 16:36:43 $

% Set the default options.
defaultopt = struct('Display','final', ...
    'Diagnostics','off', ...
    'MaxIter', '100000*numberOfVariables', ...
    'MaxRLPIter', '100*numberOfVariables', ...
    'MaxNodes', '1000*numberOfVariables', ...
    'MaxTime', 7200, ...
    'NodeDisplayInterval', 20, ...
    'TolXInteger', 1.e-8, ...
    'TolFun', 1.e-3, 'TolRLPFun', 1.e-6, ...
    'NodeSearchStrategy', 'bn', ...
    'BranchStrategy', 'maxinfeas');

% If just 'defaults' passed in, return the default options in X.
if nargin==1 && nargout <= 1 && isequal(f,'defaults')
    x = defaultopt;
    return
end

% Start timing.
startTime = cputime;
if nargin < 7
    options = [];
    if nargin < 6
        x0 = [];
        if nargin < 5
            beq = [];
            if nargin < 4
                Aeq =[];
                if (nargin < 3)
                    b = [];
                    if (nargin < 2)
                        A = [];
                    end
                end
            end
        end
    end
end

% Detect problem structure input
if nargin == 1 && isa(f,'struct')
   [f,A,b,Aeq,beq,x0,options] = separateOptimStruct(f);
end

if nargin == 0
   error(message('optim:bintprog:NotEnoughInputs'));
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(f,A,b,Aeq,beq,x0);
catch ME
    if strcmp(ME.identifier,'MATLAB:datatypes:superiorfloat')
        dataType = 'notDouble';
    end
end

if ~strcmp(dataType,'double')
    error(message('optim:bintprog:NonDoubleInput'))
end

% Options setup for the integer programming problem --------------
diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
switch optimget(options,'Display',defaultopt,'fast')
    case {'off','none'}
        verbosity = 0;
    case {'iter','iter-detailed'}
        verbosity = 2;
    case {'final','final-detailed'}
        verbosity = 1;
    otherwise
        verbosity = 1;
end

branchStrategy = optimget(options, 'BranchStrategy', defaultopt, 'fast');
nodeStrategy = optimget(options, 'NodeSearchStrategy', defaultopt, 'fast');
tolxint = optimget(options, 'TolXInteger', defaultopt, 'fast');
tolfun = optimget(options, 'TolFun', defaultopt, 'fast');
dispInterval = round( optimget(options, 'NodeDisplayInterval', defaultopt, 'fast') );
maxtime = optimget(options, 'MaxTime', defaultopt, 'fast');
maxnodes = optimget(options, 'MaxNodes', defaultopt, 'fast');
maxiter = optimget(options, 'MaxIter', defaultopt, 'fast');

% Check the input data.
[mineq, nineq] = size(A);
[meq, neq] = size(Aeq);
nvars = max([length(f) nineq neq]);

% Initialize the output argument output.
output.iterations = 0;
output.nodes = 0;
output.time = 0;
output.algorithm = 'LP-based branch-and-bound';
if strcmp(branchStrategy,'mininfeas')
  output.branchStrategy = 'minimum integer infeasibility';
else
  output.branchStrategy = 'maximum integer infeasibility';
end
if strcmp(nodeStrategy,'bn')
  output.nodeSrchStrategy = 'best node search';
else
  output.nodeSrchStrategy = 'depth first search';
end
% Check the case that the input data are empty.
if isempty(f) && isempty(A) && isempty(b) && isempty(Aeq) && isempty(beq)
    x = zeros(0,1);
    fval = 0;
    exitflag = 1;
    if verbosity >= 1
        disp(getString(message('optim:bintprog:EmptyProbEmptySol')));
    end
    return;
end

% Get the option maxiter and maxnodes after the variable nvars is defined.
if ischar(maxiter)
    if isequal(lower(maxiter),'100000*numberofvariables')
        maxiter = 100000*nvars;
    else
        error(message('optim:bintprog:InvalidFieldMaxIter'));
    end
end

if ischar(maxnodes)
    if isequal(lower(maxnodes),'1000*numberofvariables')
        maxnodes = 1000*nvars;
    else
        error(message('optim:bintprog:InvalidFieldMaxNodes'));
    end
end

% If f is empty, assume that c is a zero vector.
if isempty(f), f=zeros(nvars,1); end
if isempty(A), A=zeros(0,nvars); end
if isempty(b), b=zeros(0,1); end
if isempty(Aeq), Aeq=zeros(0,nvars); end
if isempty(beq), beq=zeros(0,1); end
if isempty(x0), x0= zeros(0,1); end

% Set vectors to column vectors
f = f(:);
b = b(:);
beq = beq(:);
x0 = x0(:);

if ~isequal(length(b),mineq)
    error(message('optim:bintprog:SizeMismatchRowsOfA'));
elseif ~isequal(length(beq),meq)
    error(message('optim:bintprog:SizeMismatchRowsOfAeq'));
elseif ~isequal(length(f),nineq) && ~isempty(A)
    error(message('optim:bintprog:SizeMismatchColsOfA'));
elseif ~isequal(length(f),neq) && ~isempty(Aeq)
    error(message('optim:bintprog:SizeMismatchColsOfAeq'));
elseif ~isequal(length(f), length(x0)) && ~isempty(x0)
    error(message('optim:bintprog:SizeMismatchLengthOfx0'));
end

if any(isnan(f)) || any(isnan(b)) || any(isnan(beq)) || any(any(isnan(A))) ...
   || any(any(isnan(Aeq)))
    error(message('optim:bintprog:InvalidInputNaN'));
end

if any(isinf(f)) || any(isinf(b)) || any(isinf(beq)) || any(any(isinf(A))) ...
   || any(any(isinf(Aeq)))
    error(message('optim:bintprog:InvalidInputInf'));
end

% Initialize the output arguments.
exitflag = -2; % Initialize exitflag as -2 (infeasible).
feasibleObjlb = -inf; % Initialize the lower bound on the objective feasibleObjlb as -inf.
feasiblePoint = false; 
feasibleX0 = false;

% Call diagnose for bintprog. 
if diagnostics
    % Do diagnostics on information so far
    gradflag = []; hessflag = []; constflag = false; gradconstflag = false; 
    non_eq=0;non_ineq=0; lin_eq=size(Aeq,1); lin_ineq=size(A,1); XOUT=ones(nvars,1);
    funfcn{1} = []; confcn{1}=[]; 
    lb = zeros(nvars,1); ub = ones(nvars, 1);
    diagnose('bintprog',output,gradflag,hessflag,constflag,gradconstflag,...
        XOUT,non_eq,non_ineq,lin_eq,lin_ineq,lb,ub,funfcn,confcn);
end

% Set up display columns for iteration on nodes. Print headers.
if verbosity >= 2
    disp(sprintf('Explored   Obj of LP   Obj of best   Unexplored   Best lower    Relative gap'));
    disp(sprintf(' nodes    relaxation  integer point     nodes    bound on obj  between bounds'));
end

% Check the feasibility of the starting point x0.
tolcon = 1.0e-6; % The default tolerance on the constraints for x0 (only).
if ~isempty(x0)
    checkineq = all(A*x0 - b <= tolcon);
    checkeq = all(abs(Aeq*x0 - beq) <= tolcon);
    checkbinary = all( abs(x0-1) < tolxint | abs(x0 -0) < tolxint );
    feasibleX0 = checkineq && checkeq && checkbinary;
    if feasibleX0
        % xinteger stores the best integer feasible point available so far.
        % If x0 is feasible, set x0 to be the initial value of xinteger and set feasibleObjective accordingly.
        xinteger = x0;
        feasibleObjective = f' * xinteger;
        feasiblePoint = feasibleX0;
        exitflag = 1; 
        % Display the given starting point x0 as nnode = 0.
        if verbosity >= 2
            disp( sprintf('*%6d   %8.4g              -     %6d             -', ...
                0, feasibleObjective, 0) );
        end
    else
        warning(message('optim:bintprog:InfeasibleStartPoint'));
    end % if feasibleX0 
end % if ~isempty(x0)

% Initialize the best integer solution if no starting point is specified
% or the given starting point is infeasible.
if ~feasiblePoint
    % xinteger stores the best integer feasible point available so far.
    % feasibleObjective is the best feasible objective value so far, initialized as inf.
    xinteger = zeros(nvars, 1);
    feasibleObjective = inf;
end
 
% Initialize a list of active nodes.
% For each node, information needed is:
% 1. nodebounds: objective of the LP relaxation of its parent node;
% 2. nodebranchvars: vector to record all branching variable indices from the root node;
% 3. nodexfix: vector to record fixed x values of -1,0,1;
% 4. nodexsol: the optimal solution of its parent node;
% 5. nodeBasicVarIdx: the optimal basis of its parent node; 
% 6. nodeNonbasicVarIdx: the nonbasic variable indices of its parent node.
if (mineq + meq > 0)
    numnodes = min(100*nvars, 10000); % Preallocation.
    if nvars <= 255
        nodebranchvars = zeros(nvars, numnodes, 'uint8');
    elseif nvars <= 65535
        nodebranchvars = zeros(nvars, numnodes, 'uint16');
    else
        error(message('optim:bintprog:BeyondRange'));
    end
    nodebounds = inf(numnodes, 1);
    nodexfix = -ones(nvars, numnodes, 'int8');  % =-1, unfixed; =0, fixed at 0; =1, fixed at 1.
    % Information of the optimal solution to the parent node.
    nodexsol = zeros(nvars, numnodes);
else % In the case that min f'*x without inequality or equality constraints. 
    nodebounds = inf;
    % Need only room for the root node.
    nodebranchvars = zeros(nvars, 1);
    nodexfix = -ones(nvars, 1);
    nodexsol = zeros(nvars, 1);
end % if (mineq + meq) > 0
nodeBasicVarIdx = []; 
nodeNonbasicVarIdx = [];

% Set up for the root node and insert it into the list.
nactivenodes = 1; % number of active nodes (generated but unexplored).
nodebounds(nactivenodes) = -inf; % Initialize the lower bound as -inf.
nodebranchvars(:, nactivenodes) = zeros(nvars,1);
nodexfix(:, nactivenodes) = -ones(nvars, 1);

% Set up options for solving LP relaxation problems by the simplex method. 
% Similar to defaultopt for linprog except that the simplex method 
% is turned on by default.
defaultoptLP = struct('Display','none',...
    'TolFun',1e-6,'Diagnostics','off',...
    'MaxIter','100*numberOfVariables', ...
    'LargeScale','off','Simplex','on');

optionsLP = defaultoptLP;
% Override LP solver default options TolFun and MaxIter
% if user specifies options.MaxRLPIter and options.TolRLPFun.
optionsLP.TolFun = optimget(options, 'TolRLPFun', defaultopt, 'fast');
optionsLP.MaxIter = optimget(options, 'MaxRLPIter', defaultopt, 'fast');
if ischar(optionsLP.MaxIter)
    if isequal(lower(optionsLP.MaxIter),'100*numberofvariables')
        optionsLP.MaxIter = 100*nvars;
    else
        error(message('optim:bintprog:InvalidValueMaxRLPIter'));
    end
end

% Set computeLambda to false as we don't need to compute Lambda in solving
% the LP relaxations.
computeLambda = false;

%-----------------------While loop ---------------
iters = 0;  % The accumulation of iterations in solving LP relaxation problems.
nexplorednodes = 0; % The number of explored nodes.
gap = inf; % gap in percentage.
while nactivenodes > 0 && gap > 100*tolfun

    % Check all the limits before exploring a new node.
    if iters >= maxiter
        exitflag = 0;   % MaxIter exceeded.
        if isinf(feasibleObjective)
            feasibleObjective = f'*xinteger;
        end
        break;
    end

    if nexplorednodes >= maxnodes
        exitflag = -4;  % MaxNodes exceeded.
        if isinf(feasibleObjective)
            feasibleObjective = f'*xinteger;
        end
        break;
    end

    if (cputime - startTime) >= maxtime
        exitflag = -5;  % MaxTime exceeded.
        if isinf(feasibleObjective)
            feasibleObjective = f'*xinteger;
        end
        break;
    end

    lb = zeros(nvars,1);
    ub = ones(nvars,1);
    % Get the next node from the active-node list as the current node.
    [feasibleObjlb, nodeindx] = getNextNode(nactivenodes, nodebounds, nodeStrategy);

    % Get the information of the current node
    xfix = nodexfix(:, nodeindx); % xfix: denotes fixed variables for the current node.
    branchvars = nodebranchvars(:, nodeindx);
    xsol = nodexsol(:, nodeindx);
    if ~isempty(nodeBasicVarIdx) 
        basicVarIdx = nodeBasicVarIdx(:, nodeindx);
        nonbasicVarIdx = nodeNonbasicVarIdx(:, nodeindx);
    else
        basicVarIdx = [];
        nonbasicVarIdx = [];
    end

    % Set up lb and ub for the current node.
    xfixindx = xfix >= 0;
    lb(xfixindx) = xfix(xfixindx);
    ub(xfixindx) = xfix(xfixindx);
    depth = length(find(branchvars));

    % Delete the node from the list of active nodes
    nodebounds(nodeindx) = [];
    nodebranchvars(:, nodeindx) = [];
    nodexfix(:, nodeindx) = [];
    nodexsol(:, nodeindx) = []; 
    if ~isempty(nodeBasicVarIdx)
        nodeBasicVarIdx(:, nodeindx) = [];
    end
    if ~isempty(nodeNonbasicVarIdx)
        nodeNonbasicVarIdx(:, nodeindx) = [];
    end
    nactivenodes = nactivenodes - 1;

    % Solve the LP relaxation of the current node.
    if ~isempty(basicVarIdx)
        % Reoptimize the LP relaxation of the current node with a dual feasible point.
        % Initialize the change of bounds on the branching variable.
        lastindxbranch = branchvars(depth);
        if xfix(lastindxbranch) == 0
            chgubval = 0;
            chgubidx = lastindxbranch;
            chglbval = [];
            chglbidx = [];
        elseif xfix(lastindxbranch) == 1
            chglbval = 1;
            chglbidx = lastindxbranch;
            chgubidx = [];
            chgubval = [];
        else
            error(message('optim:bintprog:WrongChoiceOfBranchVars'));
        end
        % Set up lb and ub for the parent node of the current node.
        lb(lastindxbranch) = 0;
        ub(lastindxbranch) = 1; 
        [cc, AA, bb, lbb, ubb, x0] = changeBoundDualFeasible(chglbidx, chglbval, chgubidx, ...
            chgubval, f, A, b, Aeq, beq, lb, ub, xsol, basicVarIdx, nonbasicVarIdx, delrows, optionsLP.TolFun);
        
        [xlp, fvallp, exitlp, outlp, basicVarIdx, nonbasicVarIdx]= dualsimplex(cc, AA, bb, lbb, ubb, x0, ...
            basicVarIdx, nonbasicVarIdx, optionsLP.MaxIter, optionsLP.TolFun, 0);
        xlp = xlp(1:nvars); % Delete slack variables.
    else % Solve the LP relaxation of the current node by simplex. 
        [xlp, fvallp, lambda, exitlp, outlp, basicVarIdx, nonbasicVarIdx, delrows] ...
            = simplex(f, A, b, Aeq, beq, lb, ub, optionsLP, defaultoptLP, computeLambda);
        if (nexplorednodes == 0) && (mineq + meq > 0) % The root node and constraints exist.
            if (nvars + mineq) <= 255
                nodeBasicVarIdx = zeros(length(basicVarIdx), numnodes , 'uint8');
                nodeNonbasicVarIdx = zeros(length(nonbasicVarIdx), numnodes, 'uint8');  
            elseif (nvars + mineq) <= 65535
                nodeBasicVarIdx = zeros(length(basicVarIdx), numnodes, 'uint16');
                nodeNonbasicVarIdx = zeros(length(nonbasicVarIdx), numnodes, 'uint16');
            else
                error(message('optim:bintprog:BeyondRange'));
            end
        end
    end % if nnz(basicVarIdx) ~= 0 
    
    nexplorednodes = nexplorednodes + 1;
    iters = iters + outlp.iterations;
    xlpisint = isint(xlp, tolxint);
    
    if exitlp == 0
        exitflag = -6; % MaxRLPIter exceeded.
        break;
    end
    
    % If the LP relaxation of the root node has an optimal solution,
    % feasibleObjlb is updated by the optimal value of it.
    if nexplorednodes == 1 && exitlp == 1
        feasibleObjlb = fvallp;
    end
     
    % Relative gap between the obj of best integer point and the best lower
    % bound on obj computed by 
    % gap (in %) = (Best Integer - Best Lower Bound)/(abs(Best Integer) + 1)*100.
    if isfinite(feasibleObjective) && isfinite(feasibleObjlb)
        gap = (feasibleObjective - feasibleObjlb)/(abs(feasibleObjective) + 1) * 100;
    end
    % Initialize the flag for new incumbent (i.e. new integer point with
    % lower objective).
    newIncumbent = false;
    
    % Pruning by the optimal value of the LP relaxation.
    % Operations depend on the solution of the current node:
    % a. if LP is infeasible or LP is optimal with greater than the obj of the best integer point,
    %    then prune the node and move to the next node;
    % b. if new incumbent (feasible integer point with lower objective value) 
    %    found, then update the current best integer node;
    % c. if LP optimal but not integer and LP optimal value less than the
    %    best integer point, then branch and add the two child nodes to the
    %    list of active nodes.
    if (exitlp == -1) || (exitlp == 1 && fvallp >= feasibleObjective)
        % Prune the node and move to the next.
        branch = false;
    elseif xlpisint && (exitlp == 1) && (fvallp < feasibleObjective) 
        % When new incumbent (i.e. integer point with lower objective)
        % arises, update the current best integer point xinteger and the
        % corresponding obj feasibleObjective.
        xinteger = xlp;
        feasibleObjective = fvallp;
        feasiblePoint = true;
        newIncumbent = true;
        exitflag = 1;
        branch = false;
        % Delete all the nodes from the active node list if their bounds >=
        % feasibleObjective because of the new incumbent.
        delindx = nodebounds(1:nactivenodes) >= feasibleObjective;   % The index of to-be-deleted nodes.
        if any(delindx)
            nodebounds(delindx) = [];
            nodebranchvars(:, delindx) = [];
            nodexfix(:, delindx) = [];
            nactivenodes = nactivenodes - nnz(delindx);
        end

        % Update the relative gap.
        gap = (feasibleObjective - feasibleObjlb)/(abs(feasibleObjective) + 1) * 100;
        
    elseif (exitlp == 1) && (fvallp < feasibleObjective) && ~xlpisint
        % Branching case: exitlp == 1 and fvallp < feasibleObjective but
        % xlp is not integer.
        fracx = xlp - floor(xlp);
        if isequal(lower(branchStrategy), 'maxinfeas')
            [t, indxbranch] = max(min(fracx, 1-fracx));
        elseif isequal(lower(branchStrategy), 'mininfeas')
            indxfrac = (fracx > tolxint);
            nonintfracx = fracx(indxfrac);
            [t, indxbranch] = max(max(nonintfracx, 1-nonintfracx));
            indx = find(indxfrac); 
            indxbranch = indx(indxbranch);
        else
            error(message('optim:bintprog:InvalidBranchStrategy'));
        end

        if xfix(indxbranch) >= 0
            error(message('optim:bintprog:WrongBranchVarSelection'));
        end

        % Branch on the left and right child nodes.
        % Insert the generated new nodes into the active node list.
        depth = depth + 1;
        branchvars(depth) = indxbranch;
   
        % Insert the right child then the left child.     
        % Set up the right child node xj = 1.
        xfix(indxbranch) = 1;
        nactivenodes = nactivenodes + 1;
        nodebounds(nactivenodes) = fvallp;
        nodebranchvars(:, nactivenodes) = branchvars;
        nodexfix(:, nactivenodes) = xfix;
        % Add the info of the optimal basis of the parent node.
        nodexsol(:, nactivenodes) = xlp;
        nodeBasicVarIdx(:,nactivenodes) = basicVarIdx;
        nodeNonbasicVarIdx(:, nactivenodes) = nonbasicVarIdx;

        % Set up the left child node.
        xfix(indxbranch) = 0;
        nactivenodes = nactivenodes + 1;
        nodebounds(nactivenodes) = fvallp;
        nodebranchvars(:, nactivenodes) = branchvars;
        nodexfix(:, nactivenodes) = xfix;
        % Add the info of the optimal basis of the parent node.
        nodexsol(:, nactivenodes) = xlp;
        nodeBasicVarIdx(:,nactivenodes) = basicVarIdx;
        nodeNonbasicVarIdx(:, nactivenodes) = nonbasicVarIdx;
    end % if (exitlp == -1) || (exitlp == 1 && fvallp >= feasibleObjective)

    % Display the iteration on nodes.
    if verbosity >= 2 && ( mod(nexplorednodes, dispInterval) == 0 || nexplorednodes == 1 || newIncumbent)
        if ~newIncumbent
            printIterOnNode(nexplorednodes, fvallp, feasibleObjective, nactivenodes, feasibleObjlb, gap);
        else  % Highlight the new incumbent.
            if isfinite(gap)
                disp( sprintf('*%6d   %10.4g   %10.4g     %6d      %10.4g    %8.2g%%', ...
                    nexplorednodes, fvallp, feasibleObjective, nactivenodes, feasibleObjlb, gap) );
            else
                disp( sprintf('*%6d   %10.4g   %10.4g     %6d      %10.4g        -', ...
                    nexplorednodes, fvallp, feasibleObjective, nactivenodes, feasibleObjlb) );
            end
        end
    end % if verbosity >= 2 && ...
    
end % while nactivenodes > 0

% Display the last iteration if it is not displayed yet.
if verbosity >= 2 && mod(nexplorednodes, dispInterval) ~= 0 && nexplorednodes ~= 1 && ~newIncumbent
   printIterOnNode(nexplorednodes, fvallp, feasibleObjective, nactivenodes, feasibleObjlb, gap); 
end

% Add a heuristic to round the result if possible.
roundx = round(xinteger);
froundx = f'*roundx;
if (froundx <= feasibleObjective) && all(A*roundx <= b) && all(Aeq*roundx == beq)
    xinteger = roundx;
    feasibleObjective = froundx;
end
    
% Throw a warning if the integer solution contains some noninteger value
% due to the integer tolerance setup.
if ~all(xinteger==1 | xinteger==0)
    warning(message('optim:bintprog:IntegerTol'));
end    

fval = feasibleObjective;
x = xinteger;
output.time = cputime - startTime;
output.nodes = nexplorednodes;
output.iterations = iters;

% Display the final statement.
if exitflag == 1
  outMessage = 'Optimization terminated.';
elseif exitflag == 0
  outMessage = 'Maximum number of iterations exceeded; increase options.MaxIter.';
elseif exitflag == -2
  outMessage = 'The problem is infeasible.';
elseif exitflag == -4
  outMessage = 'Maximum number of nodes exceeded; increase options.MaxNodes.';
elseif exitflag == -5
  outMessage = 'Maximum time exceeded; increase options.MaxTime.';
elseif exitflag == -6
  outMessage = 'Maximum number of iterations for LP Relaxation exceeded; increase options.MaxRLPIter.';
elseif exitflag ~= -1 % Note that exitflag = -1 is reserved for interrupt.
  error(message('optim:bintprog:InvalidStatus'));
end
if verbosity > 0
  disp(outMessage)
end
output.message = outMessage;

%----------------------------------------------------------------
%-----------------SUBFUNCTIONS------------------------------------
function  xisint = isint(x, tolxint)
% x is assumed to be a column vector;
% tolxint is the integrality tolerance.

if isempty(x)
    xisint = false;
else
    xisint = max(min(x-floor(x), ceil(x)- x)) <= tolxint;
end
%----------------------------------------------------------------
function printIterOnNode(nexplorednodes, fvallp, z, nactivenodes, zlb, gap)
% Print the iteration on nodes
if ~isempty(fvallp)
    if isfinite(z) && isfinite(zlb)
        disp( sprintf(' %6d   %10.4g   %10.4g     %6d      %10.4g    %8.2g%%', ...
            nexplorednodes, fvallp, z, nactivenodes, zlb, gap) );
    elseif isinf(z) && isfinite(zlb)
        disp( sprintf(' %6d   %10.4g            -     %6d      %10.4g          -', ...
            nexplorednodes, fvallp, nactivenodes, zlb) );
    elseif isfinite(z) && isinf(zlb)
        disp( sprintf(' %6d   %10.4g   %10.4g     %6d             -            -', ...
            nexplorednodes, fvallp, z, nactivenodes) );
    else
        disp( sprintf(' %6d   %10.4g           -      %6d             -            -', ...
            nexplorednodes, fvallp, nactivenodes) );
    end
else % fvallp is empty
    if isfinite(z) && isfinite(zlb)
        disp( sprintf(' %6d          -     %10.4g     %6d      %10.4g    %8.2g%%', ...
            nexplorednodes, z, nactivenodes, zlb, gap) );
    elseif isinf(z) && isfinite(zlb)
        disp( sprintf(' %6d          -              -     %6d      %10.4g          -', ...
            nexplorednodes, nactivenodes, zlb) );
    elseif isfinite(z) && isinf(zlb)
        disp( sprintf(' %6d          -     %10.4g     %6d             -            -', ...
            nexplorednodes, z, nactivenodes) );
    else
        disp( sprintf(' %6d          -              -     %6d             -            -', ...
            nexplorednodes, nactivenodes) );
    end
end
%--------------------------------------------------------------------------
function [zlb, nodeindx] = getNextNode(nactivenodes, nodebounds, nodeStrategy)
% Get the next node from the active-node list as the current node.
% Return the index of the current node and the lower bound on the objective
% function value zlb.
if isequal(nodeStrategy, 'df')
    % Choose the top node of the list as the next to be explored (nodeindx).
    nodeindx = nactivenodes;
    zlb = min(nodebounds(1:nactivenodes));
elseif isequal(nodeStrategy, 'bn')
    % Choose the node with the lowest bound from the list as the next.
    [zlb, nodeindx] = min(nodebounds(1:nactivenodes));
else
    error(message('optim:bintprog:InvalidNodeSearchStrategy'));
end
%--------------------------------------------------------------------------
%------------- END OF BINTPROG --------------------------------------------

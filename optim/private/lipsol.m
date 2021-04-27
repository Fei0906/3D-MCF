function [xsol,fval,lambda,exitflag,output] = lipsol(f,Aineq,bineq,Aeq,beq,lb,ub,options,defaultopt,computeLambda)
%LIPSOL  Linear programming Interior-Point SOLver.
%   X = LIPSOL(f,A,b) solves the linear programming problem
%
%          min f'*x    subject to:   A*x <= b
%           x
%
%   X = LIPSOL(f,A,b,Aeq,beq) solves the problem above while additionally
%   satisfying the equality constraints Aeq*x = beq.
%
%   X = LIPSOL(f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper bounds
%   on the design variables, X, so that the solution X is always in the
%   range LB <= X <= UB.  If LB is [] or if any of the entries in LB are -Inf,
%   that is taken to mean the corresponding variable is not bounded below.
%   If UB is [] or if any of the entries in UB are Inf, that is
%   taken to mean the corresponding variable is not bounded above.
%
%   X = LIPSOL(f,A,b,Aeq,beq,LB,UB,OPTIONS) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS, an 
%   argument created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolFun, LargeScale, MaxIter. Use OPTIONS = [] as a 
%   place holder if no options are set.
%
%   [X,FVAL] = LIPSOL(f,A,b) returns the value of the objective
%   function at X: FVAL = f' * X.
%
%   [X,FVAL,LAMBDA] = LIPSOL(f,A,b) returns the set of Lagrange multipliers,
%   LAMBDA, at the solution, X.
%   LAMBDA.ineqlin are the Lagrange multipliers for the inequalities.
%   LAMBDA.eqlin are the Lagrange multipliers for the equalities.
%   LAMBDA.lb the Lagrange multipliers for the lower bounds.
%   LAMBDA.ub are the Lagrange multipliers for the upper bounds.
%
%   [X,FVAL,LAMBDA,EXITFLAG] = LIPSOL(f,A,b) describes the exit conditions:
%   If EXITFLAG is:
%      > 0 then LIPSOL converged with a solution X.
%      = 0 then LIPSOL did not converge within the maximum number of iterations.
%      < 0 then the problem was infeasible.
%
%   [X,FVAL,LAMBDA,EXITFLAG,OUTPUT] = LIPSOL(f,A,b) returns a structure:
%     output.iterations: number of iterations.
%     OUTPUT.cgiterations: total number of PCG iterations.

%   Original author Yin Zhang 1995.
%   Modified by The MathWorks, Inc. since 1997.
%   Copyright 1990-2011 The University of Maryland Baltimore County and The MathWorks, Inc.
%   $Revision: 1.1.6.11 $ $Date: 2011/04/11 16:19:35 $

if (nargin < 3)
   error(message('optim:lipsol:NotEnoughInputs'))
end

n = size(f,1);
n_orig = n;
if (size(f,2) ~= 1)
   error(message('optim:lipsol:FirstInputNotCol'))
end
f_orig = f;

if (nargin < 8)
   options = [];
   if ( (nargin < 5) || isempty(Aeq) )
      Aeq = sparse(0,n);
      beq = zeros(0,1);
   end
end

if isempty(Aineq)
   Aineq = sparse(0,n);
   bineq = zeros(0,1);
end

tol = optimget(options,'TolFun',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');
switch optimget(options,'Display',defaultopt,'fast')
case {'none','off'}
   diagnostic_level = 0;
case 'iter'
   diagnostic_level = 2;
case 'final'
   diagnostic_level = 1;
case 'testing'
   diagnostic_level = 4;
otherwise
   diagnostic_level = 1;
end

if ((nargin < 7) || isempty(ub))
   ub = inf(n,1);
   if (diagnostic_level >= 4)
      fprintf('  Assuming all variables are unbounded above.\n');
   end
end
if ((nargin < 6) || isempty(lb))
   lb = -inf(n,1);
   if (diagnostic_level >= 4)
      fprintf('  Assuming all variables are unbounded below.\n');
   end
end
lbinf = ~isfinite(lb);
nlbinf = sum(lbinf);

% Collect some information for lagrange multiplier calculations
if computeLambda
   lbIn = lb;
   ubIn = ub;
   lbnotinf = ~lbinf;
   ubnotinf = isfinite(ub);
   numRowsLb = sum(lbnotinf);
   numRowsUb = sum(ubnotinf);
   numRowsAineq = size(Aineq,1);
   numRowsAeq = size(Aeq,1);
   % Indices for lagrange multipliers and initialize multipliers
   Aindex = 1:numRowsAineq+numRowsAeq;
   Aineqstart = 1; Aineqend = 0;  % end value is zero since not in A yet
   Aeqstart = 1; Aeqend = 0;      % end value is zero since not in A yet
   lbindex = 1:length(lb);
   ubindex = 1:length(ub);
   lbend = length(lb);  % we don't remove the inf ones, we just ignore them
   ubend = length(ub);  % we don't remove the inf ones, we just ignore them
   lambda.ineqlin = zeros(numRowsAineq,1);
   lambda.eqlin = zeros(numRowsAeq,1);
   lambda.upper = zeros(n_orig,1); 
   lambda.lower = zeros(n_orig,1); 
   extraAineqvars = 0;
   ilbinfubfin = [];
   Aineq_orig = Aineq; Aeq_orig = Aeq;
else
   lambda.ineqlin = [];
   lambda.eqlin = [];
   lambda.upper = [];
   lambda.lower = [];
end

% IF one of the variables is unbounded below: -Inf == lb(i) <= x(i) <= ub(i)
% THEN create a new variable splitting the old x(i) into its difference:
% x(i) -> x(i) - x(n+i) such that both are positive:
% 0 <= x(i) <= Inf and 0 <= x(n+i) <= Inf
% and a new inequality constraint reflecting the old upper bound:
% x(i) - x(n+i) <= ub(i), unless ub(i) == Inf.
% The inequality and equality constraints must be padded to reflect the new
% variable and the function f must also be adjusted.
if (nlbinf > 0)
   f = [f; -f(lbinf)];
   Aineq = [Aineq -Aineq(:,lbinf)];
   Aeq = [Aeq -Aeq(:,lbinf)];
   lbinfubfin = lbinf & isfinite(ub);
   ilbinfubfin = find(lbinfubfin);
   nlbinfubfin = sum(lbinfubfin);
   Aineq = [Aineq; ...
         sparse(1:nlbinfubfin,find(lbinfubfin),1,nlbinfubfin,n) ...
         sparse(1:nlbinfubfin,find(lbinfubfin(lbinf)),-1,nlbinfubfin,nlbinf)];
   bineq = [bineq; ub(lbinfubfin)];
   lb(lbinf) = 0;
   lb = [lb; zeros(nlbinf,1)];
   ub(lbinf) = Inf;
   ub = [ub; inf(nlbinf,1)];
   n = n + nlbinf;
   if computeLambda
      extraAineqvars = nlbinfubfin;
      Aindex = 1:length(Aindex)+ extraAineqvars;
      ubindex = 1:length(ub);
      lbindex = 1:length(lb);
   end
   if (diagnostic_level >= 4)
      fprintf(['  Splitting %d (unbounded below) variables into the' ...
            ' difference of two strictly\n  positive variables:' ...
            ' x(i) = xpos-xneg, such that xpos >= 0, xneg >= 0.\n'], nlbinf);
      if (nlbinfubfin > 0)
         fprintf(['  Adding %d inequality constraints of the form:' ...
               ' xpos - xneg <= ub(i) for the finite upper bounds.\n'], ...
            nlbinfubfin);
      end
   end
end

nslack = n;
if isempty(Aineq) % Equality constraints only; no xpos/xneg variables
   A = Aeq;
   b = beq;
   m = size(A,1);
   if computeLambda
      Aeqstart = 1; Aeqend = m;
   end
   if (diagnostic_level >= 4)
      fprintf(['  Linear Programming problem has %d equality' ...
            ' constraints on %d variables.\n'],m,n);
   end
else
   mineq = size(Aineq,1);
   meq = size(Aeq,1);
   m = mineq + meq;
   f = [f; sparse(mineq,1)];
   A = [Aineq speye(mineq); Aeq sparse(meq,mineq)];
   b = [bineq; beq];
   lb = [lb; sparse(mineq,1)];
   ub = [ub; inf(mineq,1)];
   if computeLambda
      Aineqstart = 1; Aineqend = numRowsAineq;
      Aeqstart = Aineqend + extraAineqvars + 1; Aeqend = Aeqstart + numRowsAeq - 1;
      lbindex = 1:length(lb);
      ubindex = 1:length(ub);
   end
   if (diagnostic_level >= 4)
      fprintf(['  Adding %d slack variables, one for each inequality' ...
            ' constraint, to the existing\n  set of %d variables,'],mineq,n);
   end
   n = n + mineq;
   if (diagnostic_level >= 4)
      fprintf(' resulting in %d equality constraints on %d variables.\n', ...
         mineq,n);
      if  ~isempty(Aeq)
         fprintf(['  Combining the %d (formerly inequality) constraints' ...
               ' with %d equality constraints,\n  resulting in %d equality' ...
               ' constraints on %d variables.\n'],mineq,meq,m,n);
      end
   end
end

% cholinc(A,'inf') only accepts sparse matrices.
if (~issparse(A))
   A = sparse(A);
end
f = sparse(f);

% Perform preliminary solving, rearranging and checking for infeasibility

if computeLambda
   % We might delete variables from ub/lb (and so w/z) so we need
   % a separate index into w/z 
   windex = zeros(length(ub),1);
   windex(ubindex) = 1;
   zindex = zeros(length(lb),1);
   zindex(lbindex) = 1;
   yindex = ones(size(A,1),1);
   % Save f at this point -- may be needed to compute multipliers
   f_before_deletes = f;
end

% delete fixed variables
fixed = (lb == ub);
Fixed_exist = any(fixed);
if (Fixed_exist)
   ifix = find(fixed); 
   notfixed = ~fixed;
   infx = find(notfixed);
   if computeLambda
      % Decide which bound multiplier to be nonzero according to -f(i) in case
      % we later find no other constraints are active involving this x(i);
      %  leave the other zero.
      % If we find that one of the other constraints on x(i) is active, then
      %   it is ok for both these multipliers to remain zero.
      fixedlower = ( ( f >= 0 ) & fixed );
      fixedupper = ( ( f < 0) &  fixed );
      lbindex = lbindex(infx);
      zindex = zindex(infx); 
      ubindex = ubindex(infx);
      windex = windex(infx);
   end
   xfix = lb(ifix);
   f = f(infx);
   b = b - A(:,ifix)*xfix;
   A = A(:,infx);
   lb = lb(infx);
   ub = ub(infx);
   % update n
   n = size(f,1);
   if (diagnostic_level >= 4)
       fprintf(['  Assigning %i equal values of' ...
           ' lower and upper bounds to solution.\n'],length(ifix))
   end
   % To be consistent with simplexpresolve check the right hand side, b,
   % against zero instead of a tolerance.
   if  (n == 0) && any(b ~= 0)
       theMessage = sprintf(['Exiting due to infeasibility: Violate ' ...
           'the equality constraints while all variables fixed.']);
       if (diagnostic_level >= 1)
           disp(theMessage)
       end
       exitflag = -2;
       [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
       return
   end
end

% delete zero rows
rnnzct = sum(spones(A'), 1);  
if ~isempty(A) && (any(rnnzct == 0))
   zrows = (rnnzct == 0);
   izrows = find(zrows);
   if (any(b(izrows) ~= 0))
      theMessage = sprintf(['Exiting due to infeasibility: an all-zero row in the constraint\n' ...
                     ' matrix does not have a zero in corresponding right-hand-side entry.']);
      if (diagnostic_level >= 1)
         disp(theMessage)
      end
      exitflag = -2;
      [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
      return
   else   
      nzrows = ~zrows;
      inzrows = find(nzrows);
      A = A(inzrows,:);
      b = b(inzrows); 
      rnnzct = rnnzct(inzrows);
      if (diagnostic_level >= 4)
         fprintf(['  Deleting %i all-zero rows from' ...
               ' constraint matrix.\n'],length(izrows));
      end
      if computeLambda
         % multipliers for zero rows are zero
         Aindex = Aindex(inzrows);
         yindex = yindex(inzrows);
      end
      m = size(A,1);
   end
end

% make A structurally "full rank"
sprk = sprank(A');
Row_deleted = false;
if (sprk < m) && ~isempty(A)
   Row_deleted = true;
   [dmp, tmp] = dmperm(A);
   irow = dmp(1:sprk);          % Note permutation of rows might occur 
   idelrow = dmp(sprk+1:m);
   Adel = A(idelrow,:);
   bdel = b(idelrow);
   A = A(irow,:);
   b = b(irow);
   if computeLambda
      Aindex = Aindex(irow);
      yindex = yindex(irow);
   end
   rnnzct = rnnzct(irow); 
   if (diagnostic_level >= 4)
      fprintf(['  Deleting %i dependent rows from' ...
            ' constraint matrix.\n'],m-sprk);
   end
   m = size(A,1);
end

% delete zero columns
Zrcols_exist = false;
if isempty(A)
   zrcol = 0;
else
   zrcol = (max(abs(A)) == 0)';  % max behaves differently if A only has one row:
                                 %    so zero cols won't be deleted in this case.
end
if (any(zrcol == 1))
   Zrcols_exist = true;
   izrcol = find(zrcol);
   if any( f(izrcol) < 0 & isinf(ub(izrcol)) )
      theMessage = sprintf('Exiting: the problem is unbounded.');
      if (diagnostic_level >= 1)
         disp(theMessage)
      end
      exitflag = -3; 
      [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
      return
   end
   % Set variables associated to zero columns to their optimal
   % value: either lb or ub depending on the sign of the cost 
   % coefficient. In the case where the cost coefficient is zero,
   % the variable can take any value between the bounds. We choose
   % to set it to lb  
   xzrcol = zeros(size(izrcol));   
   indx = f(izrcol)<0; 
   xzrcol(indx) = ub(izrcol(indx)); 
   indx = f(izrcol)>=0;     
   xzrcol(indx) = lb(izrcol(indx));
   
   inzcol = find(~zrcol);
   if computeLambda
      % assign a multiplier according to -f(i), leave the other zero 
      lowermultnonzero = ( f(1:n,1) >= 0 ) & zrcol(1:n,1); 
      uppermultnonzero = ( f(1:n,1) < 0 ) &  zrcol(1:n,1);  
      lambda.lower(lbindex(lowermultnonzero)) = f(lowermultnonzero); 
      lambda.upper(ubindex(uppermultnonzero)) = -f(uppermultnonzero);
      lbindex = lbindex(inzcol);                        
      ubindex = ubindex(inzcol);
      zindex = zindex(inzcol);
      windex = windex(inzcol);
   end 
   A = A(:,inzcol);
   f = f(inzcol);
   lb = lb(inzcol);                        
   ub = ub(inzcol);
   if (diagnostic_level >= 4)
      fprintf(['  Deleting %i all-zero columns from' ...
            ' constraint matrix.\n'],nnz(zrcol));
   end
   n = size(f,1);
end

% solve singleton rows
Sgtons_exist = false;
singleton = (rnnzct == 1);
nsgrows = nnz(singleton);
if (nsgrows >= max(1, .01*size(A,1)))
   Sgtons_exist = true;
   isgrows = find(singleton);
   iothers = find(~singleton);
   if (diagnostic_level >= 4)
      fprintf('  Solving %i singleton variables immediately.\n', ...
         nsgrows);
   end
   Atmp = A(isgrows,:);
   Atmp1 = spones(Atmp);
   btmp = b(isgrows);
   if (nsgrows == 1)
      isolved  = find(Atmp1); 
      insolved = find(Atmp1 == 0); 
      xsolved  = b(isgrows)/Atmp(isolved);
   else
      colnnzct = sum(Atmp1);
      isolved  = find(colnnzct);
      insolved = find(colnnzct == 0);
      [ii, jj] = find(Atmp); 
      Atmp = Atmp(ii,jj);
      btmp = btmp(ii);
      xsolved  = btmp./diag(Atmp);
      if (any(colnnzct > 1))
         repeat = diff([0; jj]) == 0;
         for i = 1:length(xsolved) - 1
            if repeat(i+1) && (xsolved(i+1) ~= xsolved(i))
               theMessage = sprintf(['Exiting due to infeasibility: singleton variables in\n' ...
                                  ' equality constraints are not feasible.']);
               if (diagnostic_level >= 1)
                  disp(theMessage)
               end
               exitflag = -2;
               [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
               return
            end
         end
         ii = find(~repeat);
         jj = ii;
         Atmp = Atmp(ii,jj);
         btmp = btmp(ii);
         xsolved  = btmp./diag(Atmp);
      end
   end
   
   % check that singleton variables are within bounds
   if any(xsolved < lb(isolved)) || any(xsolved > ub(isolved))
      theMessage = sprintf(['Exiting due to infeasibility: %i singleton variables in the equality\n' ...
                         ' constraints are not within bounds.'], ...
                        sum((xsolved<lb(isolved))|(xsolved>ub(isolved))));
      if (diagnostic_level >= 1)
         disp(theMessage)
      end
      exitflag = -2;
      [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
      return
   end
   
   if computeLambda
      % Compute which multipliers will need to be computed
      %   sgrows: what eqlin lambdas we are solving for (row in A)
      %   sgcols: what column in A that singleton is in
      sgAindex = Aindex(singleton);
      sgrows = sgAindex( (sgAindex >= Aeqstart & sgAindex <= Aeqend) ) ...
         - extraAineqvars - numRowsAineq;

      % Only want to extract out indices for the original variables, not slacks.
      %  Must also subtract out the removed variables from zero columns and fixed variables.
      [iii,jjj]=find(Atmp1');
      sgcols = lbindex( iii( iii <= (n_orig - nnz(zrcol) - nnz(fixed)) ) );
      if length(sgrows) ~= length(sgcols)
            error(message('optim:lipsol:SizeMismatch'))
      end
      Aindex = Aindex(iothers);
      yindex = yindex(iothers);
      lbindex = lbindex(insolved);
      ubindex = ubindex(insolved);
      zindex = zindex(insolved);
      windex = windex(insolved);
   end
   % reform the unsolved part of the problem
   b = b(iothers,1) - A(iothers,isolved)*xsolved;
   A = A(iothers, insolved);
   f = f(insolved);
   lb = lb(insolved);
   ub = ub(insolved);
end

n = length(ub);
nbdslack = n;
boundsInA = 0;
if (isempty(A)) % constraint matrix has been cleared out
   boundsInA = 1;
   ubfinite = find(isfinite(ub));
   lbfinite = find(isfinite(lb));
   lubfinite = length(ubfinite);
   llbfinite = length(lbfinite);
   Aineq = [sparse(1:lubfinite,ubfinite,1,lubfinite,n); ...
         sparse(1:llbfinite,lbfinite,-1,llbfinite,n)];
   % Aineqend and Aeqend are less than Aeq(Aineq)start since these are lb,ub constraints
   %   No singletons, due to slacks; no zero rows; full rank due to slacks;
   %   could have fixed variables
   bineq = [ub(ubfinite); -lb(lbfinite)];
   if (diagnostic_level >= 4)
      fprintf(['  No constraints: converting %d finite upper bounds and' ...
            ' %d finite lower bounds into %d inequality constraints.\n'], ...
         lubfinite,llbfinite,lubfinite+llbfinite);
   end
   
   mineq = size(Aineq,1);
   m = mineq;
   f = [f; sparse(mineq,1)];
   A = [Aineq speye(mineq)];
   b = bineq;
   lb = [lb; sparse(mineq,1)];
   ub = [ub; inf(mineq,1)];
   if computeLambda
      Aindex = 1:(lubfinite+llbfinite);
      Aineqstart = 1; Aineqend = 0;  % empty
      Aeqstart = 1; Aeqend = 0;  % empty
   end
   if (diagnostic_level >= 4)
      fprintf(['  Adding %d slack variables, one for each inequality' ...
            ' constraint (bounds), to the existing\n  set of %d variables,'],mineq,n);
   end
   n = n + mineq;
   if (diagnostic_level >= 4)
      fprintf(' resulting in %d equality constraints on %d variables.\n', ...
         mineq,n);
   end   
end

% shift nonzero lower bounds
Lbounds_non0 = any(lb ~= 0);
if (Lbounds_non0)
   b = b - A*lb;
   if (diagnostic_level >= 4)
      fprintf('  Shifting %i non-zero lower bounds to zero.\n', ...
         full(sum(lb ~= 0)));
   end
end
% find upper bounds
nub = 0;
iubounds = (ub ~= Inf);   
if computeLambda
   ubindex = ubindex(iubounds);  
   windex(~iubounds) = 0;        % not removing from ub/w, just zeroing out
end
Ubounds_exist = full(any(iubounds));
if (Ubounds_exist)
   ub(~iubounds) = 0;
   ub = sparse(iubounds .* (ub-lb));
   nub = nnz(ub);
end
m = size(A,1);

if computeLambda
   % now ignore the lb inf bounds
   [stmp,itmp]=setdiff(lbindex,find(lbnotinf));
   zindex(itmp) = 0;  % zero out inf bounds
   lbindex = intersect(find(lbnotinf),lbindex);
end

% Scale the problem
badscl = 1.e-4;
col_scaled = false;
absnzs = abs(nonzeros(A));
thescl = min(absnzs)/max(absnzs);

q = 1; % initialize scale factor of constraints
if (thescl < badscl)
   if (diagnostic_level >= 4)
      fprintf(['  Scaling problem by square roots of infinity norms of rows and' ...
            ' columns of constraint matrix.\n']);
   end
   
   % ----- scaling vectors ------
   absA = abs(A);
   colscl = full(sqrt(max(absA, [], 1)'));
   rowscl = full(sqrt(max(absA',[], 1)'));
   
   % ----- column scaling -----
   if (Ubounds_exist)
      ub = ub .* colscl;
   end
   
   colscl = reciprocal(colscl); 
   A = A * spdiags(colscl,0,n,n);
   
   f = f .* colscl;
   col_scaled = true;
   
   % ----- row scaling -----
   rowscl = reciprocal(rowscl);
   A = spdiags(rowscl,0,m,m) * A;
   b = b .* rowscl;
   bnrm = norm(b);
   if (bnrm > 0) 
      q = median([1 norm(f)/bnrm 1.e+8]);
      if (q > 10)
          A = q * A;
          b = q * b;
      else
          q = 1;
      end
   end
end

% Not all variables are fixed, singleton, etc, i.e. more variables left to solve for
if ~isempty(A)

    idense = [];
    ispars = 1:n;
    Dense_cols_exist = false;
    nzratio = 1;

    if (m > 2000)
        nzratio = 0.05;
    elseif (m > 1000)
        nzratio = 0.10;
    elseif (m >  500)
        nzratio = 0.20;
    end

    if (nzratio < 1)
        checking = (sum(spones(A))/m <= nzratio);
        if (any(checking == 0))
            Dense_cols_exist = true;
            % checking will be mostly ones, (1-checking) will be very sparse
            idense = find(sparse(1-checking));   % Dense  column indices
            ispars = find(checking);             % Sparse column indices
        end
    end

    if ((diagnostic_level >= 4) && ~isempty(idense))
        fprintf(['  Separating out %i dense' ...
            ' columns from constraint matrix.\n'],length(idense));
    end

    tau0 = .9995;		% step percentage to boundary
    phi0 = 1.e-5;		% initial centrality factor

    m = size(A,1);
    nt = n + nub;
    bnrm = norm(b);
    fnrm = norm(f);
    Hist = [];

    if (Ubounds_exist)
        unrm = norm(ub);
    else
        unrm = [];
    end

    if (Dense_cols_exist)
        perm = colamd(A(:,ispars)');
    else
        perm = colamd(A');
    end

    y = zeros(m,1);
    pmin = max(bnrm/100, 100);
    dmin = fnrm*.425;
    dmin = max(dmin, pmin/40);
    pmin = min(pmin, 1.e+4);
    dmin = min(dmin, 1.e+3);

    Mdense = [];
    if (Dense_cols_exist)
        Sherman_OK = true;
    else
        Sherman_OK = false;
    end

    cgiter = 0; % leaving because of backward compatibility
                % it is required for the OUTPUT.cgiterations

    % Initialize some logical parameters and factors
    ldlFactorUptodate = false;
    useNormalForm = true;

    L = [];
    Dg = [];
    Pr = [];
    Sc = [];

    % Initialize the rankDeficiency counter
    rankDeficiency = 0;
    
    % Initialize back-tracking parameters
    indStartAugForm = 1;
    lastFibonacciNum = 5;
    nextFibonacciNum = 8;

    [P,U] = getpu(A,ones(n,1),Dense_cols_exist,idense,ispars);
    if (Dense_cols_exist)
        Rinf = [];
        [temp,Sherman_OK,residual_OK,Mdense,Rinf] = ...
            ShermanSolve(P,U,full(b),1,Sherman_OK,Mdense,perm,Rinf);
        % if the Sherman-Morrison formula is not satisfactory calculate
        % with \
        if ~residual_OK
            % Disable the warnings about conditioning for singular and
            % nearly singular matrices
            warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
            warningstate2 = warning('off', 'MATLAB:singularMatrix');
            warningstate3 = warning('off', 'MATLAB:rankDeficientMatrix');
            x = A\b;
            % Restore the warning states to their original settings
            warning(warningstate1)
            warning(warningstate2)
            warning(warningstate3)
        else
            x = A' * temp;
        end
        if ~Sherman_OK
            useNormalForm = false;
        end
   else
      Rinf = builtin('_cholinf',sparse(P(perm,perm)));
      % Disable the warnings about conditioning for singular and
      % nearly singular matrices
      warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
      warningstate2 = warning('off', 'MATLAB:singularMatrix');
      temp(perm,1) = Rinf \ (Rinf' \ full(b(perm)));
      % Restore the warning states to their original settings
      warning(warningstate1)
      warning(warningstate2)

      x = A' * temp;
   end
   
   pmin = max(pmin, -min(x));
   x = max(pmin,x);
   
   z = full((f+dmin).*(f > 0) + dmin*(-dmin < f & f <= 0) - f.*(f <= -dmin));
   
   if (Ubounds_exist)
      s = spones(ub).*max(pmin, ub-x);
      w = spones(ub).*( dmin*(f > 0) + ...
         (dmin - f).*(-dmin < f & f <= 0) - 2*f.*(f <= -dmin) );
   else
      s = [];
      w = [];
   end
   
   [Rxz,Rsw,dgap] = complementy(x,z,s,w,Ubounds_exist);
   
   [Rb,Rf,rb,rf,ru,Ru] = feasibility(A,b,f,ub,Ubounds_exist,x,y,z,s,w);
   if any(isnan([rb rf ru]))
       theMessage = sprintf(['Exiting: cannot converge because the primal residual, dual residual,\n' ...
                          ' or upper-bound feasibility is NaN.']);

       if (diagnostic_level >= 1)
           disp(theMessage)
       end
       exitflag = -4; 
       [xsol, fval, output, lambda] = populateOutputs(0,0,theMessage);
       return
   end
   
   [trerror,rrb,rrf,rru,rdgap,fval,dualval] = ...
      errornobj(b,rb,bnrm,f,rf,fnrm,ub,ru,unrm,dgap,Ubounds_exist,x,y,w);
   
   iter = 0;
   converged = false;
   backed = false;
   if (~Ubounds_exist)
      sn1 = [];
   end
   
   if (diagnostic_level >= 2)
      if (Ubounds_exist)
         fprintf('\n  Residuals:   Primal     Dual     Upper    Duality     Total');
         fprintf('\n               Infeas    Infeas    Bounds     Gap        Rel\n');
         fprintf('               A*x-b   A''*y+z-w-f {x}+s-ub  x''*z+s''*w   Error\n');
         fprintf('  -------------------------------------------------------------\n');
      else
         fprintf('\n  Residuals:   Primal     Dual     Duality    Total');
         fprintf('\n               Infeas    Infeas      Gap       Rel\n');
         fprintf('               A*x-b    A''*y+z-f    x''*z      Error\n');
         fprintf('  ---------------------------------------------------\n');
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   while (iter <= maxiter)

       if (diagnostic_level >= 2)
           fprintf('  Iter %4i:  ', iter);
           if (Ubounds_exist)
               fprintf('%8.2e %8.2e %8.2e %8.2e %8.2e\n',rb,rf,ru,dgap,trerror);
           else
               fprintf('%8.2e %8.2e %8.2e %8.2e\n',rb,rf,dgap,trerror);
           end
       end

       if (iter > 0)
           stop = false;
           converged = false;
           theMessage = '';
           trerror = Hist(1,iter);
           if (trerror < tol)
               stop = true;
               converged = true;
               exitflag = 1;
               theMessage = 'Optimization terminated.';
           else
               small = 1.e-5;
               if (iter > 2)
                   blowup = (Hist(1:5,iter) > (max(tol,min(Hist(1:5,indStartAugForm:end)')')/small));
                   if (any(blowup))
                       if (useNormalForm)
                           % If a blow-up is detected when using Normal Formulation, 
                           % assign relevant quantities to their values from the 
                           % previous iteration and switch to Augmented Formulation.
                           useNormalForm = false;
                           x = xc;
                           y = yc;
                           z = zc;
                           s = sc;
                           w = wc;
                           Rxz = RxzPrevious;
                           Rsw = RswPrevious;
                           Rb = RbPrevious;
                           Rf = RfPrevious;
                           Ru = RuPrevious;
                           % Overwrite the last column of the history matrix, Hist,
                           % with the values from the iteration to which we
                           % backtracked.
                           Hist(:,end) = Hist(:,end-1);
                           indStartAugForm = iter;
                       elseif ~(useNormalForm) && (indStartAugForm==1)
                           % Reset the history if this is the first blow-up case while
                           % employing Augmented Formulation.
                           indStartAugForm = iter;
                       else
                           stop = true;
                           converged = false;
                           theMessage ...
                               = sprintf(['One or more of the residuals, duality gap, or total relative error\n', ...
                               ' has grown 100000 times greater than its minimum value so far:\n']);
                           [exitflag,msg] = detectinf(tol,Hist);
                           theMessage = [theMessage msg];
                       end % if (useNormalForm)
                   end
                   nstall = 4;
                   if (iter > nstall)
                       latest = iter-nstall:iter;
                       h = Hist(1:5,latest);
                       hmean = mean(h');
                       for i = 1:5
                           % The mean could be smaller than "small" (1e-5). Then to catch
                           % the stalling earlier also test with tol.
                           stall(i) = ((hmean(i) > min(small,tol)) & ...
                               (all(abs(h(i,:)-hmean(i)) < small*hmean(i))));
                       end
                       if (any(stall))
                           if (useNormalForm)
                               % If a stalling is detected when using Normal Formulation,
                               % assign relevant quantities to their values from a
                               % previous iteration and switch to Augmented
                               % Formulation.
                               useNormalForm = false;
                               if (iter - lastFibonacciNum) <= 6
                                   x = xHist1;
                                   y = yHist1;
                                   z = zHist1;
                                   s = sHist1;
                                   w = wHist1;
                                   Rxz = RxzHist1;
                                   Rsw = RswHist1;
                                   Rb = RbHist1;
                                   Rf = RfHist1;
                                   Ru = RuHist1;
                                   phi0 = phi0Hist1;
                                   phi = phiHist1;
                                   dgap = dgapHist1;
                               else
                                   x = xHist2;
                                   y = yHist2;
                                   z = zHist2;
                                   s = sHist2;
                                   w = wHist2;
                                   Rxz = RxzHist2;
                                   Rsw = RswHist2;
                                   Rb = RbHist2;
                                   Rf = RfHist2;
                                   Ru = RuHist2;
                                   phi0 = phi0Hist2;
                                   phi = phiHist2;
                                   dgap = dgapHist2;
                               end
                           else
                               stop = true;
                               converged = false;
                               theMessage ...
                                   = sprintf(['One or more of the residuals, duality gap, or total relative error\n', ...
                                   ' has stalled:\n']);
                               [exitflag,msg] = detectinf(tol,Hist);
                               theMessage = [theMessage msg];
                           end % if (useNormalForm)
                  end
               end
            end
         end
         if (stop)
            break % out of while (iter <= maxiter) loop
         end
      end % if (iter > 0)
      
      xn1 = reciprocal(x);
      if (Ubounds_exist)
         sn1 = reciprocal(s);
         rvmid = z.*xn1 + w.*sn1;
         vmid = reciprocal(rvmid);
      else
         rvmid = z.*xn1; 
         vmid = reciprocal(rvmid);   
      end
      
      rvmid = spdiags(rvmid,0,n,n);
      vmid = full(min(1e+15,vmid));
      
      [P,U] = getpu(A,vmid,Dense_cols_exist,idense,ispars);   
      
      [dx,dy,dz,ds,dw,ldlFactorUptodate,useNormalForm,Sherman_OK,Mdense,Rinf,L,Dg,Pr,Sc,rankDeficiency] = ...
         direction(A,P,U,Rb,Rf,Ru,Rxz,Rsw,vmid,rvmid,xn1,sn1,z,w,0,m,n,1, ...
         ldlFactorUptodate,useNormalForm,Sherman_OK,Dense_cols_exist, ...
         Ubounds_exist,Mdense,perm,Rinf,L,Dg,Pr,Sc,rankDeficiency);
      
      [ap,ad] = ratiotest(dx,dz,ds,dw,Ubounds_exist,x,z,s,w);
      
      % Saving data for backtracking
      % Blow-up
      RxzPrevious = Rxz;
      RswPrevious = Rsw;
      % Stalling
      if (iter == 1)
          RxzHist1 = Rxz;
          RswHist1 = Rsw;
      elseif (iter == nextFibonacciNum-1)
          if nextFibonacciNum > 8
              RxzHist1 = RxzHist2;
              RswHist1 = RswHist2;
          end
          RxzHist2 = Rxz;
          RswHist2 = Rsw;
      end
      
      if ((tau0*ap < 1) || (tau0*ad < 1))
         
         newdgap = (x + min(1,ap)*dx)'*(z + min(1,ad)*dz);
         if (Ubounds_exist) 
            newdgap = newdgap + (s + min(1,ap)*ds)'*(w + min(1,ad)*dw); 
         end
         sigmak = (newdgap/dgap)^2;
         sigmin = 0;
         sigmax = .208; % Zhang's choice
         p = ceil(log10(trerror)); 
         if ((p < -1) && (dgap < 1.e+3))
            sigmax = 10^(p+1);
         end
         sigmak = max(sigmin, min(sigmax, sigmak));
         mu = sigmak*dgap/nt;  
         
         Rxz = dx.*dz;
         Rsw = ds.*dw;
         
         [dx2,dy2,dz2,ds2,dw2,ldlFactorUptodate,useNormalForm,Sherman_OK,Mdense,Rinf,L,Dg,Pr,Sc,rankDeficiency] = ...
            direction(A,P,U,sparse(m,1),sparse(n,1),sparse(n,1),...
            Rxz,Rsw,vmid,rvmid,xn1,sn1,z,w,mu,m,n,2,ldlFactorUptodate,useNormalForm, ...
            Sherman_OK,Dense_cols_exist,Ubounds_exist,Mdense,perm,Rinf,L,Dg,Pr,Sc,rankDeficiency);
      
         dx = dx + dx2;
         dy = dy + dy2;
         dz = dz + dz2;
         ds = ds + ds2;
         dw = dw + dw2;
         
         [ap,ad] = ratiotest(dx,dz,ds,dw,Ubounds_exist,x,z,s,w);
         
      end
      
      tau = tau0; 
      if (~backed)
         tau = .9 + 0.1*tau0;
      end
      k = ceil(log10(trerror));
      if (k <= -5)
         tau = max(tau,1-10^k);
      end
      ap2 = min(tau*ap,1);
      ad2 = min(tau*ad,1);
      xc = x;
      yc = y;
      zc = z;
      sc = s;
      wc = w;
      % Record the variables at the first iteration and at every (Fibonacci
      % number - 1) iterations, starting from iteration 7.
      if (iter == 1)
          xHist1 = x;
          yHist1 = y;
          zHist1 = z;
          sHist1 = s;
          wHist1 = w;
          phi0Hist1 = phi0;
          phiHist1 = phi;
          dgapHist1 = dgap;
      elseif (iter == nextFibonacciNum-1)
          if nextFibonacciNum > 8
              xHist1 = xHist2;
              yHist1 = yHist2;
              zHist1 = zHist2;
              sHist1 = sHist2;
              wHist1 = wHist2;
              phi0Hist1 = phi0Hist2;
              phiHist1 = phiHist2;
              dgapHist1 = dgapHist2;
          end
          xHist2 = x;
          yHist2 = y;
          zHist2 = z;
          sHist2 = s;
          wHist2 = w;
          phi0Hist2 = phi0;
          phiHist2 = phi;
          dgapHist2 = dgap;
      end
      step = [1 .9975 .95 .90 .75 .50];
      for k = 1:length(step)
         x = xc + step(k)*ap2*dx;
         y = yc + step(k)*ad2*dy;
         z = zc + step(k)*ad2*dz;
         s = sc + step(k)*ap2*ds;
         w = wc + step(k)*ad2*dw;
         [Rxz,Rsw,dgap] = complementy(x,z,s,w,Ubounds_exist);
         phi = nt*full(min([Rxz;Rsw(find(Rsw))]))/dgap;
         if ((max(ap2,ad2) == 1) || (phi >= phi0))
            break
         end
      end
      phi0 = min(phi0, phi);
      if (k > 1)
          backed = true;
      end
      % Saving data for backtracking and moving up in the Fibonacci
      % sequence
      % Blow-up
      RbPrevious = Rb;
      RfPrevious = Rf;
      RuPrevious = Ru;
      % Stalling
      if (iter == 1)
          RbHist1 = Rb;
          RfHist1 = Rf;
          RuHist1 = Ru;
      elseif (iter == nextFibonacciNum-1)
          if nextFibonacciNum > 8
              RbHist1 = RbHist2;
              RfHist1 = RfHist2;
              RuHist1 = RuHist2;
          end
          RbHist2 = Rb;
          RfHist2 = Rf;
          RuHist2 = Ru;
          tmpNum = nextFibonacciNum;
          nextFibonacciNum = lastFibonacciNum + nextFibonacciNum;
          lastFibonacciNum = tmpNum;
      end
      
      [Rb,Rf,rb,rf,ru,Ru] = feasibility(A,b,f,ub,Ubounds_exist,x,y,z,s,w);
      if any(isnan([rb rf ru]))
          theMessage = sprintf(['Exiting: cannot converge because the primal residual, dual residual,\n' ...
                         ' or upper-bound feasibility is NaN.']);
          if (diagnostic_level >= 1)
              disp(theMessage)
          end
          exitflag = -4; 
          [xsol, fval, output, lambda] = populateOutputs(iter,cgiter,theMessage);
          return
      end
      
      [trerror,rrb,rrf,rru,rdgap,fval,dualval] = ...
         errornobj(b,rb,bnrm,f,rf,fnrm,ub,ru,unrm,dgap,Ubounds_exist,x,y,w);
      
      Hist = [Hist [trerror rrb rrf rru rdgap fval dualval]'];
      iter = iter + 1;
      
   end % while (iter <= maxiter)
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
else % All variables are fixed, singleton, etc - no more variables left to solve for.
   iter = 0;
   cgiter = 0;
   exitflag = 1;
   theMessage = 'Solution determined by constraints.';    
   converged = true;
   x = zeros(0,1);
   y = zeros(0,1);
   z = zeros(0,1);
   w = zeros(0,1);
   s = [];
end

xsol = x;

if (col_scaled)
   xsol = xsol .* colscl;
   if (diagnostic_level >= 4)
      fprintf('  Scaling solution back to original problem.\n');
   end
end

if (Lbounds_non0)
   xsol = xsol + lb;
   if (diagnostic_level >= 4)
      fprintf(['  Shifting solution back away from zero to reflect' ...
            ' original lower bounds.\n']);
   end
end

if (nbdslack < n)
   % Remove slack solution variables
   xsol = xsol(1:nbdslack);
   if (diagnostic_level >= 4)
      fprintf('  Eliminating %d slack variables (bounds) from solution.\n', ...
         n-nbdslack);
   end
end

if (Sgtons_exist)
   xtmp(insolved,1) = xsol;
   xtmp(isolved,1) = xsolved;
   xsol = xtmp; 
   n = length(xsol);
   if (diagnostic_level >= 4)
      fprintf('  Reuniting %d singleton variables with solution.\n', ...
         length(xsolved));
   end
end

if (Zrcols_exist)
   xtmp(inzcol,1) = xsol;
   xtmp(izrcol,1) = xzrcol;
   xsol = xtmp;
   n = length(xsol);
   if (diagnostic_level >= 4)
      fprintf('  Reuniting %d obvious zero variables with solution.\n', ...
         length(xzrcol));
   end
end

if (Row_deleted)
   if (norm(Adel*xsol-bdel)/max(1,norm(bdel)) > tol)
      converged = false;
      exitflag = -2;
      theMessage = sprintf(['The primal is infeasible; the equality constraints are dependent\n' ...
                         ' but not consistent.']);
   else
      if (diagnostic_level >= 4)
         fprintf('  Solution satisfies %d dependent constraints.\n', ...
            size(Adel,1));
      end
   end
end

if (Fixed_exist)
   xtmp(infx,1) = xsol;
   xtmp(ifix,1) = xfix;
   xsol = xtmp;
   n = length(xsol);
   if (diagnostic_level >= 4)
      fprintf('  Reuniting %d fixed variables with solution.\n', ...
         length(ifix));
   end
end

if (nslack < n)
   % Remove slack solution variables
   xsol = xsol(1:nslack);
   if (diagnostic_level >= 4)
      fprintf('  Eliminating %d slack variables from solution.\n', ...
         n-nslack);
   end
end

if (nlbinf > 0) % if there were unbounded below variables split up
   xsol(lbinf) = xsol(lbinf) - xsol(n_orig+1:end);
   xsol = xsol(1:n_orig);
   if (diagnostic_level >= 4)
      fprintf(['  Subtracting the strictly positive differences of %d' ...
            ' unbounded below variables.\n'], ...
         nlbinf);
   end
end

fval = full(f_orig' * xsol);

if computeLambda
   % Calculate lagrange multipliers
   lbindexfinal = lbindex(lbindex <= lbend); 
   ubindexfinal = ubindex(ubindex <= ubend);
   zindexfinal = find(zindex);
   windexfinal = find(windex);
   nnzlbindexfinal = nnz(lbindexfinal);
   nnzubindexfinal = nnz(ubindexfinal);
   if boundsInA
      lambda.upper(ubindexfinal) = w(windexfinal(1:nnzubindexfinal)) - y(1:nnzubindexfinal);
      lambda.lower(lbindexfinal) = z(zindexfinal(1:nnzlbindexfinal)) - ...
         y(nnzubindexfinal+1:nnzlbindexfinal+nnzubindexfinal);  
      if Fixed_exist 
         lambda.lower(fixedlower) = f_before_deletes(fixedlower);
         lambda.upper(fixedupper) = -f_before_deletes(fixedupper);
      end
      if Sgtons_exist && ~isempty(sgrows)
         % sgrows: what eqlin lambdas we are solving for (row in A)
         % sgcols: what column in A that singleton is in
         result = -f_orig + lambda.lower - lambda.upper;
         lambda.eqlin(sgrows) = (result(sgcols)) ./ diag(Aeq_orig(sgrows, sgcols));
      end
   else
      iAineqindex = Aindex( find (Aindex >= Aineqstart & Aindex <= Aineqend) );
      % Subtract off extraAineqvars and numRowsAineq even if those were
      %    removed since the indices still reflect their existence
      iAeqindex = Aindex(find (Aindex >= Aeqstart & Aindex <= Aeqend) ) ...
         - extraAineqvars - numRowsAineq;
      
      nnzAineqindex = nnz(iAineqindex);
      nnzAeqindex = nnz(iAeqindex);
      
      % In case scaling, transform the dual variables to the original system 
      if (col_scaled)  
          y = y .* rowscl * q; 
          if ~isempty(z)
              z =  z ./ colscl;
          end
          if ~isempty(w)
              w =  w ./ colscl;
          end
      end
      
      lambda.ineqlin(iAineqindex) = full(-y(1:nnzAineqindex));
      yAeq = y(nnzAineqindex+extraAineqvars+1:end);
      lambda.eqlin(iAeqindex) = full(-yAeq(1:nnzAeqindex));
      
      nnzubindex = nnz(ubindexfinal);
      lambda.upper(ubindexfinal) = full(w(windexfinal(1:nnzubindexfinal)));
      % if we recast some upper bound constraints as xpos/xneg constraints
      extraindex = Aindex((Aineqend+1):(Aeqstart-1)); 
      lambda.upper(ilbinfubfin) = full(-y(extraindex));
      lambda.lower(lbindexfinal) = full(z(zindexfinal(1:nnzlbindexfinal)));
      
      if Fixed_exist && ~Sgtons_exist
         Alambda = Aineq_orig'*lambda.ineqlin + Aeq_orig'*lambda.eqlin;
         lambda.lower(fixedlower) = f_orig(fixedlower) + Alambda(fixedlower);
         lambda.upper(fixedupper) = - f_orig(fixedupper) - Alambda(fixedupper);
         
      elseif Sgtons_exist && ~Fixed_exist && ~isempty(sgrows)
         % sgrows: what eqlin lambdas we are solving for (row in A)
         % sgcols: what column in A that singleton is in
         result = -f_orig - Aineq_orig'*lambda.ineqlin - Aeq_orig'*lambda.eqlin ...
            + lambda.lower - lambda.upper;
         lambda.eqlin(sgrows) = (result(sgcols))./ diag(Aeq_orig(sgrows, sgcols));
      elseif Fixed_exist && Sgtons_exist  
         % solve for singletons first
         if ~isempty(sgrows)
             result = -f_orig - Aineq_orig'*lambda.ineqlin - Aeq_orig'*lambda.eqlin ...
                 + lambda.lower - lambda.upper;
             lambda.eqlin(sgrows) = (result(sgcols))./ diag(Aeq_orig(sgrows, sgcols));
         end
         % now the fixed variables multipliers
         Alambda = Aineq_orig'*lambda.ineqlin + Aeq_orig'*lambda.eqlin;
         lambda.lower(fixedlower) = f_orig(fixedlower) + Alambda(fixedlower);
         lambda.upper(fixedupper) = - f_orig(fixedupper) - Alambda(fixedupper);
      end
   end
   % Keep multipliers corresponding to the original variable bounds only.
   lambda.lower = lambda.lower(1:n_orig);
   lambda.upper = lambda.upper(1:n_orig);
end  % ComputeLambda

% Counter "iter" incremented at very end of main loop; when maxiter
% reached, "iter" is one unit higher that number of iterations performed.
if ~converged && iter == maxiter+1
   theMessage = sprintf('Maximum number of iterations exceeded; increase options.MaxIter.');
   exitflag = 0;
   iter = iter - 1; 
end   

if ~converged
  theMessage = sprintf('Exiting: %s', theMessage');
end

if (diagnostic_level >= 1)
  disp(theMessage)
end

output.iterations = iter;
output.algorithm = 'large-scale: interior point';
output.cgiterations = cgiter;
output.message = theMessage;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIPSOL Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exitflag,msg] = detectinf(tol,Hist)
%DETECTINF  Used to detect an infeasible problem.
%   One or more of these 5 convergence criteria either blew up or stalled:
%   trerror  total relative error      max([rrb rrf rru rdgap])
%   rrb      relative primal residual  norm(A*x-b) / max(1,norm(b))
%   rrf      relative dual residual    norm(A'*y+z-w-f) / max(1,norm(f))
%   rru      relative upper bounds     norm(spones(s).*x+s-ub) / max(1,norm(ub))
%   rdgap    relative objective gap    abs(fval-dualval) / max([1 abs(fval) abs(dualval)])
%   where
%   fval     primal objective          f'*x
%   dualval  dual objective            b'*y - w'*ub      
%
%   Note: these convergence criteria are stored in the first 5 rows of Hist,
%   while fval and dualval are the 6th and 7th rows of Hist respectively.

% minimum of all relative primal residuals (take upper bounds into account)
minrp = min(Hist(2,:)+Hist(4,:));
% minimum of all relative dual residuals
minrd = min(Hist(3,:)); 

if minrp < tol
   exitflag = -3;
   msg = sprintf(['         the dual appears to be infeasible (and the primal unbounded).' ,...
              '      \n         (The primal residual < TolFun=%3.2e.)'],tol);
   return
end
if minrd < tol
   exitflag = -2;
   msg = sprintf(['         the primal appears to be infeasible (and the dual unbounded).', ...
              '\n         (The dual residual < TolFun=%3.2e.)'],tol);
   return
end

tol10 = 10 * tol;
sqrttol = sqrt(tol);
if ((minrp < tol10) && (minrd > sqrttol))
   exitflag = -3;
   msg = sprintf(['         the dual appears to be infeasible (and the primal unbounded) since', ...
         '\n         the dual residual > sqrt(TolFun)=%3.2e.' ,...
         '\n         (The primal residual < 10*TolFun=%3.2e.)'], ...
      sqrttol,tol10);
   return
end
if ((minrd < tol10) && (minrp > sqrttol))
   exitflag = -2;
   msg = sprintf(['         the primal appears to be infeasible (and the dual unbounded) since', ...
             '\n         the primal residual > sqrt(TolFun)=%3.2e.', ...
             '\n         (The dual residual < 10*TolFun=%3.2e.)'], ...
          sqrttol,tol10);
   return
end

iter = size(Hist,2);
fval = Hist(6,iter);
dualval = Hist(7,iter);
if ((fval < -1.e+10) && (dualval < 1.e+6))
    exitflag = -3;
    msg = sprintf(['         the dual appears to be infeasible and the primal unbounded since' ,...
          '\n         the primal objective < -1e+10',...
          '\n         and the dual objective < 1e+6.']);
   return
end
if ((dualval > 1.e+10) && (fval > -1.e+6))
   exitflag = -2;
   msg = sprintf(['         the primal appears to be infeasible and the dual unbounded since' ,...
                   '\n         the dual objective > 1e+10',...
          '\n         and the primal objective > -1e+6.']);
   return
end

exitflag = -5;
msg = sprintf('         both the primal and the dual appear to be infeasible.\n');
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rxz,Rsw,dgap] = complementy(x,z,s,w,Ubounds_exist)
%COMPLEMENTY  Evaluate the complementarity vectors and gap.

Rxz = x .* z;
if (Ubounds_exist)
   Rsw = s .* w;
else
   Rsw = [];
end
dgap = full(sum([Rxz; Rsw]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rb,Rf,rb,rf,ru,Ru] = feasibility(A,b,f,ub,Ubounds_exist,x,y,z,s,w)
%FEASIBILITY  Evaluate feasibility residual vectors and their norms.

Rb = A*x - b;
Rf = A'*y + z - f;
if (Ubounds_exist)
   Rf = Rf - w;
   Ru = spones(s).*x + s - ub;
else
   Ru = [];
end
rb = norm(Rb);
rf = norm(Rf);
ru = norm(Ru);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [trerror,rrb,rrf,rru,rdgap,fval,dualval] = ...
   errornobj(b,rb,bnrm,f,rf,fnrm,ub,ru,unrm,dgap,Ubounds_exist,x,y,w)
%ERRORNOBJ  Calculate the total relative error and objective values.

rrb = rb/max(1,bnrm);
rrf = rf/max(1,fnrm);
rru = 0;
fval = full(f'*x);
dualval = full(b'*y);
if (Ubounds_exist) 
   rru = ru/(1+unrm); 
   dualval = dualval - w'*ub;
end
if (min(abs(fval),abs(dualval)) < 1.e-4)
   rdgap = abs(fval-dualval);
else
   rdgap = abs(fval-dualval)/max([1 abs(fval) abs(dualval)]);
end
trerror = max([rrb rrf rru rdgap]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = reciprocal(X)
%RECIPROCAL  Invert the nonzero entries of a matrix elementwise.
%   Y = RECIPROCAL(X) has the same sparsity pattern as X
%	 (except possibly for underflow).

if (issparse(X))
   [m, n]  = size(X);
   [i,j,Y] = find(X);
   Y = sparse(i,j,1./Y,m,n);
else
   Y = 1./X;
end

Y = min(Y,1e12); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,U] = getpu(A,vmid,Dense_cols_exist,idense,ispars)
%GETPU  Compute Matrix P and U

[m,n] = size(A);

if (Dense_cols_exist)
   ns = length(ispars);
   nd = length(idense); 
   Ds = spdiags(vmid(ispars),0,ns,ns);
   P = A(:,ispars) * Ds * A(:,ispars)';
   Dd = spdiags(sqrt(vmid(idense)),0,nd,nd);
   U = full(A(:,idense) * Dd);
else
   P = A * spdiags(vmid,0,n,n) * A';
   U = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dx,dy,dz,ds,dw,ldlFactorUptodate,useNormalForm,Sherman_OK,Mdense,Rinf,L,Dg,Pr,Sc,rankDeficiency] = ...
   direction(A,P,U,Rb,Rf,Ru,Rxz,Rsw,vmid,rvmid,xn1,sn1,z,w,mu,m,n,flag,...
   ldlFactorUptodate,useNormalForm,Sherman_OK,Dense_cols_exist, ...
   Ubounds_exist,Mdense,perm,Rinf,L,Dg,Pr,Sc,rankDeficiency)
%DIRECTION  Compute search directions

if (mu ~= 0)
   Rxz = Rxz - mu;
end
Rf = Rf - Rxz .* xn1;
if (Ubounds_exist)
   if (mu ~= 0)
      Rsw = Rsw - mu;
   end
   Rf = Rf + (Rsw - Ru .* w) .* sn1;
end

if (useNormalForm)
    rhs = -(Rb + A * (vmid .* Rf));
    ldlFactorUptodate = false;
    if (~Dense_cols_exist)
        if (flag == 1)
            Rinf = builtin('_cholinf',sparse(P(perm,perm)));
        end
        % Disable the warnings about conditioning for singular and
        % nearly singular matrices
        warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
        warningstate2 = warning('off', 'MATLAB:singularMatrix');
        dy(perm,1) = Rinf \ (Rinf' \ full(rhs(perm)));
        % Restore the warning states to their original settings
        warning(warningstate1)
        warning(warningstate2)
        dx = vmid .* (A' * dy + Rf);
    elseif (Sherman_OK) % Dense_cols_exist and Sherman_OK
        [dy,Sherman_OK,residual_OK,Mdense,Rinf] = ...
            ShermanSolve(P,U,full(rhs),flag,Sherman_OK,Mdense,perm,Rinf);
        if (Sherman_OK) % residual_OK is also true
            dx = vmid .* (A' * dy + Rf);
        else % if ~Sherman_OK
            useNormalForm = false;
            if (residual_OK)
                dx = vmid .* (A' * dy + Rf);
                dz = -(z .* dx + Rxz) .* xn1;
                if (Ubounds_exist)
                    ds = -(dx .* spones(w) + Ru);
                    dw = -(w .* ds + Rsw) .* sn1;
                else
                    ds = [];
                    dw = [];
                end
                return 
            %else continue and enter the next if-clause
            end % if (residual_OK)
        end % if (Sherman_OK)
    % Else, Dense_cols_exist and ~Sherman_OK,
    % can only happen outside the main loop, because whenever
    % ~Sherman_OK in this function useNormalForm is assigned to false and
    % the if-clause with useNormalForm as the condition is never entered.
    %
    end % if (~Dense_cols_exist)
end % if (useNormalForm)

if ~(useNormalForm) % use Augmented Formulation instead.
    if (flag == 1) || ((flag == 2) && ~ldlFactorUptodate)
        % If the pedictor iteration or the corrector iteration without
        % up-to-date ldl factors, perform the factorization.
        [L,Dg,Pr,Sc,~,matRank] = ldl([-rvmid A';A sparse(m,m)],'vector'); 
        if (matRank < m+n)
            rankDeficiency = m+n-matRank;
            Dg(end-rankDeficiency+1:end,end-rankDeficiency+1:end) = 1.0;
        else
            rankDeficiency = 0;
        end        
        ldlFactorUptodate = true;
    end
    rhs = [-Rf; -Rb];
    zeroassignVector = ones(m+n,1);
    zeroassignVector(end-rankDeficiency+1:end) = 0.0; 
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
    warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
    warningstate2 = warning('off', 'MATLAB:singularMatrix');  
    dxdy = Sc(:,Pr)*(L'\(Dg\(L\(zeroassignVector.*(Sc(Pr,:)*rhs)))));    
    % Restore the warning states to their original settings
    warning(warningstate1)
    warning(warningstate2)
    dx = dxdy(1:n);
    dy = dxdy(n+1:end);
end % if ~(useNormalForm)

dz = -(z .* dx + Rxz) .* xn1;
if (Ubounds_exist)
   ds = -(dx .* spones(w) + Ru);
   dw = -(w .* ds + Rsw) .* sn1;
else
   ds = [];
   dw = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ap,ad] = ratiotest(dx,dz,ds,dw,Ubounds_exist,x,z,s,w)
%RATIOTEST  Ratio test

ap = -1/min([dx(find(x))./nonzeros(x); -0.1]);
ad = -1/min([dz(find(z))./nonzeros(z); -0.1]);
if (Ubounds_exist)
   as = -1/min([ds(find(s))./nonzeros(s); -0.1]);
   aw = -1/min([dw(find(w))./nonzeros(w); -0.1]);
   ap = min(ap, as);
   ad = min(ad, aw);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,Sherman_OK,residual_OK,Mdense,Rinf] = ...
    ShermanSolve(P,U,b,flag,Sherman_OK,Mdense,perm,Rinf)
%ShermanSolve  Solve linear system with dense columns using the
%    Sherman-Morrison formula and check the solution residuals.

bnrm = norm(b);
x = zeros(size(b));
tol1  = min(1.e-2, 1.e-7*(1+bnrm));
tol2  = 10*tol1;
residual_OK = false;

if (Sherman_OK)
   [x,Mdense,Rinf] = sherman(P,U,b,flag,Mdense,perm,Rinf);
end
resid = norm(b - (P*x + U*(U'*x)));
if (resid > tol1)
   Sherman_OK = false;
end
if (resid < tol2)
   residual_OK = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,Mdense,Rinf] = sherman(P,U,b,flag,Mdense,perm,Rinf)
%SHERMAN  Use the Sherman-Morrison formula to "solve" (P + U*U')x = b
%   where P is sparse and U is low-rank.

% Disable the warnings about conditioning for singular and
% nearly singular matrices
warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
warningstate2 = warning('off', 'MATLAB:singularMatrix');

if (flag == 1)
   nd = size(U,2); 
   PinvU = zeros(size(U));
   Rinf = builtin('_cholinf',sparse(P(perm,perm)));
   PinvU(perm,:) = Rinf \ (Rinf' \ U(perm,:));
   Mdense = eye(nd) + U' * PinvU;
end

x(perm,1) = Rinf \ (Rinf' \ b(perm));
tmp = U * (Mdense \ (U' * x));
tmp(perm,1) = Rinf \ (Rinf' \ tmp(perm));
% Restore the warning states to their original settings
warning(warningstate1)
warning(warningstate2)
x = x - tmp;

%-----------------------------------------------------------------------------
function [xsol, fval, output, lambda] = populateOutputs(iter,cgiter,theMessage)
%POPULATEOUTPUTS sets all of the output values (except exitflag)
% to the standard values in an early termination of LINPROG.

xsol = [];
fval = [];

lambda.ineqlin = [];
lambda.eqlin = [];
lambda.upper = [];
lambda.lower = [];

output.iterations = iter;
output.algorithm = 'large-scale: interior point';
output.cgiterations = cgiter;
output.message = theMessage;
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% End of LIPSOL and its subfunctions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

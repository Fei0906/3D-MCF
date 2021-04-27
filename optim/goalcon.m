function [GS,GSeq,GGS,GGSeq] = goalcon(V,neqgoal,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin)
%

%GOALCON Utility function to translate gradient in goal-attainment problem.
%   Intermediate function used to translate goal attainment
%   problem into constrained optimization problem.
%   Used by FGOALATTAIN and FMINIMAX.
%
%   See also GOALFUN.

%   Copyright 1990-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.11 $  $Date: 2011/06/30 16:36:48 $

numVar = length(V) - 1;
x(:) = V(1:numVar);
lambda = V(numVar + 1);

% Compute the constraints
f = []; gf = [];  % Tell parser that f and g are variables.
gc = []; gceq = []; c = []; ceq = [];

% Determine whether or not to check for errors in objective functions
if ~errCheck         
    switch funfcn{1} % evaluate objective functions and possibly gradients
        case 'fun'
            f = feval(funfcn{3},x,varargin{:});  
        case 'fungrad'
            [f,gf] = feval(funfcn{3},x,varargin{:});
        case 'fun_then_grad'
            f = feval(funfcn{3},x,varargin{:});  
            gf = feval(funfcn{4},x,varargin{:});
        otherwise
            error(message('optim:goalcon:UndefinedCalltypeFgoalattain'))
    end

    % Evaluate constraints
    switch confcn{1}
        case 'fun'
            [ctmp,ceqtmp] = feval(confcn{3},x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
            gc = [];
            gceq = [];
        case 'fungrad'
            [ctmp,ceqtmp,gc,gceq] = feval(confcn{3},x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
        case 'fun_then_grad'
            [ctmp,ceqtmp] = feval(confcn{3},x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
            [gc,gceq] = feval(confcn{4},x,varargin{:});
        case ''
            c = []; ceq = [];
            gc = [];
            gceq = [];
        otherwise
            error(message('optim:goalcon:UndefinedCalltypeFgoalattain'))
    end

    nGoalCnstr = numel(f); % number of goal attainment constraints 
else
    % Evaluate objective functions and check for internal errors
    switch funfcn{1} % evaluate function and possibly gradients
        case 'fun'
            try
                f = feval(funfcn{3},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:ObjectiveError', ...
                    getString(message('optim:goalcon:ObjectiveError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
        case 'fungrad'
            try
                [f,gf] = feval(funfcn{3},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:ObjectiveError', ...
                    getString(message('optim:goalcon:ObjectiveError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
        case 'fun_then_grad'
            try
                f = feval(funfcn{3},x,varargin{:});              
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:ObjectiveError', ...
                    getString(message('optim:goalcon:ObjectiveError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            try
                gf = feval(funfcn{4},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:GradientError', ...
                    getString(message('optim:goalcon:GradientError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
        otherwise
            error(message('optim:goalcon:UndefinedCalltypeFgoalattain'));
    end

    % Evaluate constraints and check for internal errors
    switch confcn{1}
        case 'fun'
            try
                [ctmp,ceqtmp] = feval(confcn{3},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:NonlconError', ...
                    getString(message('optim:goalcon:NonlconError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            c = ctmp(:); ceq = ceqtmp(:);
            gc = [];
            gceq = [];
        case 'fungrad'
            try
                [ctmp,ceqtmp,gc,gceq] = feval(confcn{3},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:NonlconError', ...
                    getString(message('optim:goalcon:NonlconError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            c = ctmp(:); ceq = ceqtmp(:);
        case 'fun_then_grad'
            try
                [ctmp,ceqtmp] = feval(confcn{3},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:NonlconError', ...
                    getString(message('optim:goalcon:NonlconError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
            c = ctmp(:); ceq = ceqtmp(:);
            try
                [gc,gceq] = feval(confcn{4},x,varargin{:});
            catch userFcn_ME
                optim_ME = MException('optim:goalcon:NonlconGradError', ...
                    getString(message('optim:goalcon:NonlconGradError')));
                userFcn_ME = addCause(userFcn_ME,optim_ME);
                rethrow(userFcn_ME)
            end
        case ''
            c = []; ceq = [];
            gc = [];
            gceq = [];
        otherwise
            error(message('optim:goalcon:UndefinedCalltypeFgoalattain'));
    end

    nGoalCnstr = numel(f); % number of goal attainment constraints 

    % Initial check of user-functions. The errCheck flag is set accordingly
    % in calling functions. Check the size of user-supplied objective
    % function gradients
    if ~isempty(gf) && ~all(size(gf) == [numVar,nGoalCnstr])
        gradSizeExcept = MException('optim:goalcon:InvalidSizeOfGF', ...
            getString(message('optim:goalcon:InvalidSizeOfGF',numVar,nGoalCnstr)));
        throwAsCaller(gradSizeExcept)
    end

    % Check that the length of user objective function is equal to length(GOAL)
    if length(GOAL) ~= nGoalCnstr
        goalSizeExcept = MException('optim:goalcon:InvalidGoalAndFunSizes', ...
            getString(message('optim:goalcon:InvalidGoalAndFunSizes')));
        throwAsCaller(goalSizeExcept)
    end
    
    non_eq = length(ceq);
    non_ineq = length(c);

    [cgrow, cgcol]= size(gc);
    [ceqgrow, ceqgcol]= size(gceq);
    
    % Check the size of user-supplied non-linear constraint gradients
    if ~isempty(gc) && (cgrow ~= numVar || cgcol ~= non_ineq)
        cGradSizeExcept = MException('optim:goalcon:InvalidSizeOfGC', ...
            getString(message('optim:goalcon:InvalidSizeOfGC',numVar,non_ineq)));
        throwAsCaller(cGradSizeExcept)
    end
    if ~isempty(gceq) && (ceqgrow ~= numVar || ceqgcol ~= non_eq)
        ceqGradSizeExcept = MException('optim:goalcon:InvalidSizeOfGCeq', ... 
            getString(message('optim:goalcon:InvalidSizeOfGCeq',numVar,non_eq)));
        throwAsCaller(ceqGradSizeExcept)
    end
end

% Calculate goal attainment constraints ( f(i)-GOAL(i) )/WEIGHT(i)-lambda <= 0
GS = zeros(nGoalCnstr+neqgoal,1);
for i=1:nGoalCnstr
     if WEIGHT(i)~=0
       diff=f(i)-GOAL(i);
       GS(i)=sign(real(diff))*norm(diff)/WEIGHT(i)-lambda; 
       if i<=neqgoal % neqgoal comes from options.GoalsExactAchieve 
          GS(i+nGoalCnstr)=-GS(i)-2*lambda; % f(i)+lambda*WEIGHT>=GOAL
       end
     else % hard constraint
       GS(i)=f(i)-GOAL(i);
       if i<=neqgoal 
          GS(i+nGoalCnstr)=-GS(i)-1e-10; % f(i)>=GOAL(i)-1e-10
       end
     end
end

% Append goal attainment constraint at the end of inequalities vector
GS=[c(:);GS]; 

% Equalities and gradient matrix of equalities
GSeq = ceq(:);
size_ceq = size(GSeq);
GGSeq = [gceq; zeros(1,size_ceq(1))]; % add zero row (derivatives w.r.t. lambda)

if isempty(gf) && isempty(gc)
   GGS=[]; GGSeq = [];
elseif (isempty(gf) && ~isempty(gc)) 
   error(message('optim:goalcon:GradFunRequiredForGradCon'));
else % grad of f available, grad of inequalities available or not
   
   % Gradient matrix of inequality constraints: 
   % [grad_x of inequalities c(x)          | grad_x of attainment constraint]
   % [zero row (derivatives w.r.t. lambda) | -1 -1 . . . . . .  . . . . ..-1]
   GL = -ones(1,nGoalCnstr+neqgoal);
   for i=1:nGoalCnstr
      if WEIGHT(i)~=0
         gf(:,i)=gf(:,i)/WEIGHT(i);
         if i<=neqgoal,
            gf(:,i+nGoalCnstr)=-gf(:,i);
         end 
      else % hard constraint
         GL(1,i)=0; 
      end
   end
   
   GGS=[gf;GL];
   
   sizegc=size(gc);
   % Put gc first
   if sizegc(1)>0, 
      GGS=[[gc;zeros(1,sizegc(2))],GGS]; 
   end
end

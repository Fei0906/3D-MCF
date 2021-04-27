function [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s,numEvals] = finDiffFseminf(xCurrent,...
             xOriginalShape,funfcn,lb,ub,fCurrent,cCurrent,options, ...
             variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,startnlineq,varargin)
%

%finDiffFseminf computes finite-difference derivatives for fseminf.
%
% This helper function computes finite-difference derivatives of the objective 
% and constraint functions.
%
%  [gradf,cJac] = finDiffFseminf(...) 
%
% computes the finite-difference gradients of the objective and
% constraint functions.
%
%  gradf = finDiffFseminf(...)
%
% computes the finite-difference gradients of the objective function.
%
% INPUT:
% xCurrent              Point where gradient is desired
% xOriginalShape        Shape of the vector of variables supplied by the user
%                       (The value of xOriginalShape is NOT used)
% funfcn, confcn        Cell arrays containing info about objective and
%                       constraints, respectively. The objective (constraint) 
%                       derivatives are computed if and only if funfcn 
%                       (confcn) is nonempty.
%                       
% lb, ub                Lower and upper bounds
% fCurrent, cCurrent    Values at xCurrent of the function and the constraints 
%                       to be differentiated. Note that fCurrent can be a scalar 
%                       or a (row or column) vector. 
%
% options               Structure containing option for the finite difference calculation.
%  Fields:
%  DiffMinChange,        
%  DiffMaxChange        Minimum and maximum values of perturbation of xCurrent 
%  FinDiffType          'forward' or 'central' for forward or central finite
%                       differences, respectively
%  TypicalX             Magnitude of a typical value of variable x. 
%  FinDiffRelStep       Perturbation factor to compute the perturbation
%
% variables             Variables w.r.t which we want to differentiate. Possible 
%                       values are 'all' or an integer between 1 and the
%                       number of variables.
%
% LAMBDA,NEWLAMBDA,
% OLDLAMBDA,POINT,
% FLAG,s,startnlineq    Parameters for semi-infinite constraints
%
% varargin              Problem-dependent parameters passed to the objective and 
%                       constraint functions
%
% OUTPUT:
% gradf                 If fCurrent is a scalar, gradf is the finite-difference 
%                       gradient of the objective; if fCurrent is a vector,
%                       gradf is the finite-difference Jacobian  
% cJac                  Finite-difference Jacobian of the constraints
% NEWLAMBDA,
% OLDLAMBDA,s           Parameters for semi-infinite constraints. 

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:05:48 $

DiffMinChange = options.DiffMinChange;
DiffMaxChange = options.DiffMaxChange;
FinDiffRelStep = options.FinDiffRelStep;
typicalx = options.TypicalX;
finDiffType = options.FinDiffType;
% fwdFinDiff is true for 'forward' and false for 'central'. Default to 'forward'.
fwdFinDiff = strcmpi(finDiffType,'forward') || isempty(finDiffType);

% For vector-valued functions in funfcn, we make the function values 
% (both at xCurrent and those at the perturbed points) column vectors 
% to ensure that the given fCurrent and the computed fplus will have the 
% same shape so that fplus - fCurrent be well defined.  
fCurrent = fCurrent(:);
numberOfFunctions = numel(fCurrent);
numberOfVariables = numel(xCurrent); 
functionIsScalar = (numberOfFunctions == 1);

% nonEmptyLowerBounds = true if lb is not empty, false if it's empty;
% analogoulsy for nonEmptyUpperBound
nonEmptyLowerBounds = ~isempty(lb);
nonEmptyUpperBounds = ~isempty(ub);

% Make sure xCurrent and typicalx are column vectors so that the 
% operation max(abs(xCurrent),abs(typicalx)) won't error
xCurrent = xCurrent(:); typicalx = typicalx(:);
if fwdFinDiff % forward differences
    CHG = FinDiffRelStep.*nonzerosign(xCurrent).*max(abs(xCurrent),abs(typicalx));
else % central differences
    CHG = FinDiffRelStep.*max(abs(xCurrent),abs(typicalx));
end
%
% Make sure step size lies within DiffminChange and DiffMaxChange
%
CHG = nonzerosign(CHG).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
len_cCurrent = length(cCurrent); % for semi-infinite

if nargout < 3
   NEWLAMBDA=[]; OLDLAMBDA=[]; s=[];
end
if nargout > 1
      cJac = zeros(len_cCurrent,numberOfVariables);  
else
      cJac = [];
end
% allVariables = true/false if finite-differencing wrt to all/one variables
allVariables = false;
if ischar(variables)
   if strcmp(variables,'all')
      variables = 1:numberOfVariables;
      allVariables = true;
   else
      error(message('optim:fseminf:InvalidVariables'))
   end
end

% Preallocate gradf for speed 
if ~isempty(funfcn)
    if functionIsScalar
        gradf = zeros(numberOfVariables,1);
    elseif allVariables % vector-function and gradf estimates full Jacobian
        gradf = zeros(numberOfFunctions,numberOfVariables);
    else % vector-function and gradf estimates one column of Jacobian
        gradf = zeros(numberOfFunctions,1);
    end
else
    gradf = [];
end

% Do this switch outside of loop for speed
vararginAdjusted = varargin(3:end);

if fwdFinDiff % forward differences
    for gcnt = variables
        if (nonEmptyLowerBounds && isfinite(lb(gcnt))) || (nonEmptyUpperBounds && isfinite(ub(gcnt)))
            % Enforce bounds while finite-differencing; may need to modify step
            CHG(gcnt) = fwdFinDiffInsideBnds(xCurrent(gcnt),lb(gcnt),ub(gcnt),CHG(gcnt),gcnt,DiffMinChange);
        end 
        
        % Save and perturb component of current point
        xCElementOrig = xCurrent(gcnt);
        xCurrent(gcnt) = xCElementOrig + CHG(gcnt);
        
        xOriginalShape(:) = xCurrent;
        if ~isempty(funfcn) % objective gradient required
            fplus = feval(funfcn{3},xOriginalShape,vararginAdjusted{:});
        end
        
        if ~isempty(funfcn)
            if functionIsScalar
                gradf(gcnt,1) =  (fplus(:)-fCurrent)/CHG(gcnt);
            elseif allVariables % vector-function and gradf estimates full Jacobian
                gradf(:,gcnt) = (fplus(:)-fCurrent)/CHG(gcnt);
            else % vector-function and gradf estimates only one column of Jacobian
                gradf = (fplus(:)-fCurrent)/CHG(gcnt);
            end
        end
        
        if ~isempty(cJac) % constraint gradient required
            if gcnt == numberOfVariables,
                FLAG = -1;
            end
            [ctmp,ceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,s] = ...
                semicon(xOriginalShape,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,startnlineq,varargin{:});
            cplus = [ceqtmp(:); ctmp(:)];
            
            % Next line used for problems with varying number of constraints
            if len_cCurrent~=length(cplus)
                cplus = v2sort(cCurrent,cplus);
            end
            if ~isempty(cplus)
                cJac(:,gcnt) = (cplus - cCurrent)/CHG(gcnt);
            end
        end
        xCurrent(gcnt) = xCElementOrig; % restore perturbed component of current point
    end % for
    numEvals = numel(variables); % number of evaluations
else % central differences
    for gcnt = variables
        if (nonEmptyLowerBounds && isfinite(lb(gcnt))) || (nonEmptyUpperBounds && isfinite(ub(gcnt)))
            % Enforce bounds while finite-differencing; may need to modify step
            [CHG(gcnt),formulaType] = cntrlFinDiffInsideBnds(xCurrent(gcnt),lb(gcnt),ub(gcnt),CHG(gcnt),gcnt,DiffMinChange);
        else
            formulaType = 0; % regular central differences if no bounds
        end

        xCElementOrig = xCurrent(gcnt); % save component of current point
        % Perturb component of current point
        if formulaType == 0      % central difference formula
            xPert1Element = xCElementOrig - CHG(gcnt);
            xPert2Element = xCElementOrig + CHG(gcnt);
        elseif formulaType == -1 % double left step difference formula
            xPert1Element = xCElementOrig - 2*CHG(gcnt);
            xPert2Element = xCElementOrig - CHG(gcnt);
        else                     % double right step difference formula
            xPert1Element = xCElementOrig + CHG(gcnt);
            xPert2Element = xCElementOrig + 2*CHG(gcnt);
        end
        
        % Form first perturbed point
        xCurrent(gcnt) = xPert1Element;
        xOriginalShape(:) = xCurrent;
        % Evaluate objective at first perturbed point
        if ~isempty(funfcn) % objective gradient required
            fPert1 = feval(funfcn{3},xOriginalShape,vararginAdjusted{:});
        end
        % Evaluate nonlinear constraints at first perturbed point
        if ~isempty(cJac) % constraint gradient required
            if gcnt == numberOfVariables,
                FLAG = -1;
            end
            [ctmp,ceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,s] = ...
                semicon(xOriginalShape,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,startnlineq,varargin{:});
            cPert1 = [ceqtmp(:); ctmp(:)];
        end
        % Next line used for problems with varying number of constraints
        if len_cCurrent~=length(cPert1)
            cPert1 = v2sort(cCurrent,cPert1);
        end
        
        % Form second perturbed point
        xCurrent(gcnt) = xPert2Element;
        xOriginalShape(:) = xCurrent;
        % Evaluate objective at second perturbed point
        if ~isempty(funfcn) % objective gradient required
            fPert2 = feval(funfcn{3},xOriginalShape,vararginAdjusted{:});
        end
        % Evaluate nonlinear constraints at second perturbed point
        if ~isempty(cJac) % constraint gradient required
            if gcnt == numberOfVariables,
                FLAG = -1;
            end
            [ctmp,ceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,s] = ...
                semicon(xOriginalShape,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,startnlineq,varargin{:});
            cPert2 = [ceqtmp(:); ctmp(:)];
        end
        % Next line used for problems with varying number of constraints
        if len_cCurrent~=length(cPert2)
            cPert2 = v2sort(cCurrent,cPert2);
        end
        
        % Compute finite difference gradients
        if ~isempty(funfcn)
            gradf(gcnt,1) = twoStepFinDiffFormulas(formulaType,CHG(gcnt),fCurrent,fPert1,fPert2);
        end
        if ~isempty(cJac)
            if ~isempty(cPert1) && ~isempty(cPert2)
                cJac(:,gcnt) = twoStepFinDiffFormulas(formulaType,CHG(gcnt),cCurrent,cPert1,cPert2);
            end
        end
        
        xCurrent(gcnt) = xCElementOrig; % restore perturbed component of current point
    end
    numEvals = 2*numel(variables); % number of evaluations
end








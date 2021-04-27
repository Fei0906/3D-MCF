function [step,quadObj,normStep,normStepScal] = ...
    dogleg(nvar,F,JAC,grad,Delta,scalMat,varargin)
%DOGLEG approximately solves trust region subproblem via a dogleg approach.
%
%   DOGLEG finds an approximate solution d to the problem:
%
%     min_d      f + g'd + 0.5*d'Bd
%
%     subject to ||Dd|| <= Delta
%
%   where g is the gradient of f, B is a Hessian approximation of f, D is a
%   diagonal scaling matrix and Delta is a given trust region radius.

%   Copyright 1990-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2010/12/03 23:18:46 $

% NOTE: The scaling matrix D above is called scalMat in this routine.

% Compute scaled gradient and other scaled terms.
gradscal = grad./scalMat;
gradscal2 = gradscal./scalMat;
normgradscal = norm(gradscal);

if normgradscal >= eps
    % First compute the Cauchy step (in scaled space).
    dCauchy = -(Delta/normgradscal)*gradscal;
    JACvec = JAC*gradscal2;
    denom = Delta*(JACvec'*JACvec);
    tauterm = normgradscal^3/denom;
    tauC = min(1,tauterm);
    dCauchy = tauC*dCauchy;
    
    % Compute quadratic objective at Cauchy point.
    JACvec = JAC*(dCauchy./scalMat); 
    objCauchy = gradscal'*dCauchy + 0.5*(JACvec'*JACvec);   
    normdCauchy = min(norm(dCauchy),Delta);
else
    % Set Cauchy step to zero step and continue.
    objCauchy = 0;
    normdCauchy = 0;
    dCauchy = zeros(nvar,1);
end

if Delta - normdCauchy < eps;
    % Take the Cauchy step if it is at the boundary of the trust region.
    step = dCauchy; quadObj = objCauchy;
else
    % Compute the Gauss-Newton step (in scaled space).
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
    warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
    warningstate2 = warning('off', 'MATLAB:singularMatrix');
    dNewton = -JAC\F;
    % Restore the warning states to their original settings
    warning(warningstate1)
    warning(warningstate2)
    
    dNewton = dNewton.*scalMat;     % scale the step
    
    if any(~isfinite(dNewton))
        % Take the Cauchy step if the Gauss-Newton step gives bad values.
        step = dCauchy; quadObj = objCauchy;
    else 
        normdNewt = norm(dNewton);
        if normdNewt <= Delta
            % Use the Newton direction as the trial step
            step = dNewton;
        else
            % Find the intersect point along dogleg path.
            Delta2 = Delta^2;
            normdCauchy2 = min(normdCauchy^2,Delta2);
            normdNewt2 = normdNewt^2;
            dCdN = dCauchy'*dNewton;
            dCdNdist2 = max((normdCauchy2+normdNewt2-2*dCdN),0);
            
            if dCdNdist2 == 0
                tauI = 0;
            else
                tauI = (normdCauchy2-dCdN + sqrt((dCdN-normdCauchy2)^2 ...
                    + dCdNdist2*(Delta2-normdCauchy2))) / dCdNdist2;
            end
            step = dCauchy + tauI*(dNewton-dCauchy);
        end
        % Compute quadratic objective at trial point.
        JACvec = JAC*(step./scalMat);
        quadObj = gradscal'*step + 0.5*(JACvec'*JACvec);
        
        % Compare Cauchy step and trial step (Newton or intersection)
        if objCauchy < quadObj
            step = dCauchy; quadObj = objCauchy;
        end
    end
end

% The step computed was the scaled step.  Unscale it.
normStepScal = norm(step);
step = step./scalMat;
normStep = norm(step);





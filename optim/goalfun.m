function [lambda,gf]=goalfun(V,neqcstr,FUNfcns,GRADfcns,WEIGHT,GOAL,x,sizeCheck,varargin)
%

%GOALFUN Utility function to translate goal-attainment problem.
%
%   Intermediate function used to translate goal attainment
%   problem into constrained optimization problem.
%   Used by FGOALATTAIN, FMINIMAX and FMINCON (from private/NLCONST). 
%   
%   See also GOALCON.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/09 01:05:56 $

nx = length(V) - 1;

% Compute the objective function 
lambda=V(nx + 1);   

if nargout > 1
   % Compute the gradient of the objective
   gf=[zeros(nx,1); 1];
end

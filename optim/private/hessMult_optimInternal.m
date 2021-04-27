function W = hessMult_optimInternal(Hinfo,Y,varargin)
%hessMult_optimInternal	Hessian-matrix product
%
% An example of a Hessian-matrix product function
% where Hinfo is the actual Hessian and so W = Hinfo*Y.
%
% Note: varargin is not used but must be provided in case 
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/05/07 18:25:12 $

W = Hinfo*Y;

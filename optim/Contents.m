% Optimization Toolbox
% Version 6.2 (R2012a) 29-Dec-2011 
%
% Nonlinear minimization of functions.
%   fminbnd      - Scalar bounded nonlinear function minimization.
%   fmincon      - Multidimensional constrained nonlinear minimization.
%   fminsearch   - Multidimensional unconstrained nonlinear minimization, 
%                  by Nelder-Mead direct search method.
%   fminunc      - Multidimensional unconstrained nonlinear minimization.
%   fseminf      - Multidimensional constrained minimization, semi-infinite 
%                  constraints.
%   ktrlink      - Multidimensional constrained nonlinear minimization 
%                  using KNITRO(R) third-party libraries.
%
% Nonlinear minimization of multi-objective functions.
%   fgoalattain  - Multidimensional goal attainment optimization 
%   fminimax     - Multidimensional minimax optimization.
%        
% Linear least squares (of matrix problems).
%   lsqlin       - Linear least squares with linear constraints.
%   lsqnonneg    - Linear least squares with nonnegativity constraints.
%
% Nonlinear least squares (of functions).
%   lsqcurvefit  - Nonlinear curvefitting via least squares (with bounds).
%   lsqnonlin    - Nonlinear least squares with upper and lower bounds.
%
% Nonlinear zero finding (equation solving).
%   fzero        - Scalar nonlinear zero finding.
%   fsolve       - Nonlinear system of equations solve (function solve).
%
% Minimization of matrix problems.
%   bintprog     - Binary integer (linear) programming.
%   linprog      - Linear programming.
%   quadprog     - Quadratic programming.
%
% Controlling defaults and options.
%   optimset     - Create or alter optimization OPTIONS structure. 
%   optimget     - Get optimization parameters from OPTIONS structure. 
%
% Utility Routines.
%   color        - Column partition for sparse finite differences.
%   fzmult       - Multiplication with fundamental nullspace basis.
%   gangstr      - Zero out 'small' entries subject to structural rank.
%
% Graphical user interface and plot routines
%   optimtool                   - Optimization Toolbox Graphical User 
%                                 Interface
%   optimplotconstrviolation    - Plot max. constraint violation at each 
%                                 iteration
%   optimplotfirstorderopt      - Plot first-order optimality at each 
%                                 iteration
%   optimplotresnorm            - Plot value of the norm of residuals at
%                                 each iteration
%   optimplotstepsize           - Plot step size at each iteration

%   Copyright 1990-2011 The MathWorks, Inc.
%   Generated from Contents.m_template revision 1.1.6.5  $Date: 2009/04/15 23:21:22 $

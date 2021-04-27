function stop = optimplotresnorm(x,optimValues,state,varargin)
% OPTIMPLOTRESNORM Plot value of the norm of residuals at each iteration.
%
%   STOP = OPTIMPLOTRESNORM(X,OPTIMVALUES,STATE) plots OPTIMVALUES.resnorm.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTRESNORM as the plot
%   function
%     options = optimset('PlotFcns',@optimplotresnorm);
%
%   Pass the options into an optimization problem to view the plot
%     lsqnonlin(@(x) sin(3*x),[1 4],[],[],options);

%   Copyright 2006-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/06/30 16:36:54 $

persistent plotavailable
stop = false;

switch state
     case 'init'
         if isfield(optimValues,'resnorm')
             plotavailable = true;
         else
             plotavailable = false;
             title(getString(message('optim:optimplot:TitleNormResid', ...
                    getString(message('optim:optimplot:NotAvailable')))),'interp','none');
         end
    case 'iter'
        if plotavailable
            if optimValues.iteration == 0
                % The 'iter' case is  called during the zeroth iteration,
                % but it now has values that were empty during the 'init' case
                plotresnorm = plot(optimValues.iteration,optimValues.resnorm,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                xlabel(getString(message('optim:optimplot:XlabelIter')),'interp','none');
                set(plotresnorm,'Tag','optimplotresnorm');
                ylabel(getString(message('optim:optimplot:YlabelNormResid')),'interp','none');
                title(getString(message('optim:optimplot:TitleNormResid', ...
                    sprintf('%g',norm(optimValues.resnorm)))),'interp','none');
            else
                plotresnorm = findobj(get(gca,'Children'),'Tag','optimplotresnorm');
                newX = [get(plotresnorm,'Xdata') optimValues.iteration];
                newY = [get(plotresnorm,'Ydata') optimValues.resnorm];
                set(plotresnorm,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String', ...
                    getString(message('optim:optimplot:TitleNormResid', ...
                    sprintf('%g',optimValues.resnorm))));
            end
        end
end

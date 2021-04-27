function stop = optimplotfirstorderopt(x,optimValues,state,varargin)
% OPTIMPLOTFIRSTORDEROPT Plot first-order optimality at each iteration.
%
%   STOP = OPTIMPLOTFIRSTORDEROPT(X,OPTIMVALUES,STATE) plots
%   OPTIMVALUES.firstorderopt.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTFIRSTORDEROPT as the
%   plot function
%     options = optimset('PlotFcns',@optimplotfirstorderopt);
%
%   Pass the options into an optimization problem to view the plot
%      fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0],[],[],options)

%   Copyright 2006-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2011/06/30 16:36:53 $

persistent plotavailable
stop = false;

switch state
    case 'iter'
        if optimValues.iteration == 1
            if isfield(optimValues,'firstorderopt') && ~isempty(optimValues.firstorderopt)
                plotavailable = true;

                % The 'iter' case is  called during the zeroth iteration, but
                % firstorderopt may still  be empty.  Start plotting at the
                % first iteration.
                plotfirstorderopt = plot(optimValues.iteration,optimValues.firstorderopt,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                title(getString(message('optim:optimplot:TitleFirstOrderOpt', ...
                    sprintf('%g',optimValues.firstorderopt))),'interp','none');
                xlabel(getString(message('optim:optimplot:XlabelIter')),'interp','none');
                ylabel(getString(message('optim:optimplot:YlabelFirstOrderOpt')),'interp','none');
                set(plotfirstorderopt,'Tag','optimplotfirstorderopt');
            else % firstorderopt field does not exist or is empty
                plotavailable = false;
                title(getString(message('optim:optimplot:TitleFirstOrderOpt', ...
                    getString(message('optim:optimplot:NotAvailable')))),'interp','none');
            end
        else
            if plotavailable
                plotfirstorderopt = findobj(get(gca,'Children'),'Tag','optimplotfirstorderopt');
                newX = [get(plotfirstorderopt,'Xdata') optimValues.iteration];
                newY = [get(plotfirstorderopt,'Ydata') optimValues.firstorderopt];
                set(plotfirstorderopt,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String', ...
                    getString(message('optim:optimplot:TitleFirstOrderOpt', ...
                    sprintf('%g',optimValues.firstorderopt))));
            end
        end
end
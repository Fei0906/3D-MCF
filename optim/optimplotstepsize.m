function stop = optimplotstepsize(x,optimValues,state,varargin)
% OPTIMPLOTSTEPSIZE Plot step size at each iteration.
%
%   STOP = OPTIMPLOTSTEPSIZE(X,OPTIMVALUES,STATE) plots 
%   OPTIMVALUES.stepsize.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTSTEPSIZE as the plot
%   function
%     options = optimset('PlotFcns',@optimplotstepsize);
%
%   Pass the options into an optimization problem to view the plot
%      fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0],[],[],options)

%   Copyright 2006-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2011/06/30 16:36:55 $

persistent plotavailable
stop = false;

switch state
    case 'init'
        if isfield(optimValues,'stepsize')
            plotavailable = true;
        else
            plotavailable = false;
            title(getString(message('optim:optimplot:TitleStepSize', ...
                getString(message('optim:optimplot:NotAvailable')))),'interp','none');
        end
    case 'iter'
        if plotavailable
            if optimValues.iteration == 1
                % The 'iter' case is  called during the zeroth iteration, but
                % stepsize is still empty.  Start plotting at the first
                % iteration.
                plotstepsize = plot(optimValues.iteration,optimValues.stepsize,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                title(getString(message('optim:optimplot:TitleStepSize', ...
                    sprintf('%g',optimValues.stepsize))),'interp','none');
                xlabel(getString(message('optim:optimplot:XlabelIter')),'interp','none');
                ylabel(getString(message('optim:optimplot:YlabelStepSize')),'interp','none');
                set(plotstepsize,'Tag','optimplotstepsize');
            else
                plotstepsize = findobj(get(gca,'Children'),'Tag','optimplotstepsize');
                newX = [get(plotstepsize,'Xdata') optimValues.iteration];
                newY = [get(plotstepsize,'Ydata') optimValues.stepsize];
                set(plotstepsize,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String', ...
                    getString(message('optim:optimplot:TitleStepSize', ...
                    sprintf('%g',optimValues.stepsize))));
            end
        end
end

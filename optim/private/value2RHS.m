function RHS = value2RHS(value)
%value2RHS Convert a value into a valid MATLAB expression.

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2009/08/08 01:11:59 $



if(length(value) == 1)
    if(isa(value,'function_handle'))
        RHS = func2str(value);
        % Anonymous function add @ by default.
        if RHS(1) ~= '@'
          RHS = ['@' RHS];
        end
    elseif(isnumeric(value))
        RHS = num2str(value);
        % If value is double with no decimal digits, e.g., 1.0 then num2str 
        % converts to char and not a string. taking care of this case explicitly
        if isscalar(RHS)
           RHS = [RHS ' '];
        end
    elseif(ischar(value))
        RHS = ['''' value ''''];
    elseif(iscell(value))
        RHS = sprintf('{ %s }',value2RHS(value{1}));
    elseif(isa(value,'inline'))
        RHS = sprintf('inline(''%s'')',char(value));
    elseif (isa(value,'struct'))
        RHS = '<userStructure>';
    elseif isobject(value)
        RHS = '<userClass>';
    else
        RHS = '<userData>';
    end
elseif(length(value) > 1)
    if(isnumeric(value))
        RHS = '[ ';
        r = size(value,1);
        for i = 1:r
            RHS = [RHS, num2strNoSpace(value(i,:))]; % we use num2strNoSpace because num2str puts more spaces than we need
            if i~=r
                RHS = [RHS , ' ; '];
            end    
        end
        RHS = [ RHS ' ]' ];
    elseif(iscell(value))
        RHS = '{ ';
        [r, c] = size(value);
        for i = 1:r
            for j = 1:c
                RHS = [RHS, ' ', value2RHS(value{i,j}) ];
            end
            if i ~= r
                RHS = [RHS , ' ; '];
            end
        end
        RHS = [ RHS ' }' ];
    elseif ischar(value) && size(value,1) > 1% it must be an array of strings
         RHS = '[ ';
        r = size(value,1);
        for i = 1:r
                RHS = [RHS, '''' value(i,:) ''''];
            if i ~= r
                RHS = [RHS , ' ; '];
            end
        end
        RHS = [ RHS ' ]' ];
    elseif ischar(value) && size(value,1) == 1% A string!
        RHS = ['''' value ''''];
    elseif (isa(value,'struct'))
        RHS = '<userStructure>';
    elseif isobject(value)
        RHS = '<userClass>';
    else
        RHS = '<userData>';
    end
else % it's empty!
    RHS = '[]';
end
%------------------------------------------------------------- 
function stringWithoutSpaces = num2strNoSpace(value)
% Convert to char array
stringWithSpaces = num2str(value);
% Remove extra white spaces
stringWithoutSpaces = '';
while (true)
   [token, stringWithSpaces] = strtok(stringWithSpaces);
   stringWithoutSpaces = [stringWithoutSpaces,' ', token];
   if isempty(stringWithSpaces)
       return;
   end
end

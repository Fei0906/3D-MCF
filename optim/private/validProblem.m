function [isValid,errmsg] = validProblem(myStruct,requiredFields,validValues)
%validProblem validates optimization problem structure.
%   Checks if 'myStruct' contains all the fields in the cell
%   array 'requiredFields'. To validate the values of 'requiredFields'
%   pass 'validValues'. The argument 'validValues' must be a nested cell
%   array of valid values. The output argument 'isValid' is a  boolean.
%
%   Example:
%
%    % Create a problem structure that has wrong value for 'solver' field
%      probStruct = struct('solver','fminx','options',optimset);
%    % Suppose requiredFields are 'solver' and 'options'
%    % Get valid solver names using createProblemStruct and assume options
%    % have no known valid values (any value is okay)
%      validValues = {fieldnames(createProblemStruct('solvers')), {} };
%    % Validate the structure 'probStruct'
%      [isValid,errmsg] = validProblem(probStruct,{'solver','options'},validValues)

%   Private to OPTIMTOOL

%   Copyright 2005-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2011/06/30 16:37:14 $


errmsg = '';
for i = 1:length(requiredFields)
    field = requiredFields{i};
    % Check the field name is valid
    if ~ismember(field,fieldnames(myStruct))
        errmsg = [errmsg, getString(message('optim:optimtool:ProbStructMissingField',field))];
        continue;
    end
    % Now check the values of the field
    okayValues = validValues{i};
    if isempty(okayValues) % no values to compare; valid values
        continue;
    else
        validValue = false;
        for j = 1:length(okayValues)  % check against valid values
            if isequal(myStruct.(field),okayValues{j})
                validValue = true;
                break;
            end
        end
        if ~validValue 
            errmsg = [errmsg, getString(message('optim:optimtool:InvalidFieldInProbStruct',field))];
        end
    end
end
isValid = isempty(errmsg);
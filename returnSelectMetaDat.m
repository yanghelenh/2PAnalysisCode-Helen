% returnSelectMetaDat.m
%
% Function that takes in metaDat struct, metaDat variable, and condition to
%  filter on for the metaDat variable and returns indices of valid data as
%  well as metaDat struct with only the valid rows.
%
% Operates on either single variable or cell array of them. For numerical
% variables, condition to filter on should be single relational operation
% that can be evaluated with eval.
%
% Function doesn't error check for correct syntax in filtering
%
% INPUTS:
%   metaDat - metadata struct, output of loadMetadataSpreadsheet.m
%   vars - metaDat variable(s) to filter on
%   conds - conditions to filter on
%
% OUTPUTS:
%   inds - indicies of selected data (into metaDat struct's fields
%   filtMetaDat - metaDat struct with only selected fields
%
% CREATED: 8/24/19 - HHY
%
% UPDATED: 8/26/19 - HHY
%

function [inds, filtMetaDat] = returnSelectMetaDat(metaDat, vars, conds)

    % preallocate logical array, start with all elements of metaDat
    logicalDat = ones(size(metaDat.ExperimentName));

    % multiple conditions
    if iscell(vars)
        % loop through all conditions
        for i = 1:length(vars)
            % deal with age comparison
            if strcmpi(vars{i}, 'Age')
                logOut = zeros(size(metaDat.Age));
                for j = 1:length(metaDat.Age)
                    output = eval([metaDat.Age{j} conds{i}]);
                    logOut(j) = logical(prod(output));
                end
            % if the field is one of the numerical fields
            elseif strcmpi(class(metaDat.(vars{i})), 'double')
                logOut = eval([metaDat.(vars{i}) conds{i}]);
            % field is one of the string fields
            else
                logOut = strcmpi(metaDat.(vars{i}), conds{i});
            end
            
            % combine with previous 
            logicalDat = logicalDat .* logOut;
        end
    else % only one condition
        if strcmpi(vars, 'Age')
            logicalDat = zeros(size(metaDat.Age));
            for j = 1:length(metaDat.Age)
                output = eval([metaDat.Age{j} conds]);
                logicalDat(j) = logical(prod(output));
            end
        elseif strcmpi(class(metaDat.(vars)), 'double')
            logicalDat = eval([metaDat.(vars) conds]);
        else
            logicalDat = strcmpi(metaDat.(vars), conds);
        end
    end
    
    % convert logical to indicies
    inds = find(logicalDat);
    
    % metaDat with only selected data
    metaDatFields = fieldnames(metaDat);
    
    % loop through all fields, copy them into filtMetaDat, only selected
    %  data
    for i = 1:length(metaDatFields)
        filtMetaDat.(metaDatFields{i}) = ...
            metaDat.(metaDatFields{i})(logicalDat);
    end

end
% plotVarExpl.m
%
% 


function plotVarExpl(kernels)
    % get all fields in kernels struct (which kernels computed)
    kernelFields = fieldnames(kernels);
    numFields = length(kernelFields);
    numROIs = length(kernels);
    
    numElm = length(kernelFields)*length(kernels);
    
    ind = 1:numElm;
    
    
    c = 1;
    allVarExpl = zeros(size(ind));
    for i = 1:numFields
        for j = 1:numROIs
            allVarExpl(c) = kernels(j).(kernelFields{i}).varExpl;
            fprintf('%d: %s %d \n', c, kernelFields{i}, j);
            c = c+1;
        end
    end
    
    figure;
    scatter(ind,allVarExpl, 'filled');
    ylim([0 0.5]);
    
    
end
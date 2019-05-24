% plotAllKernelPred.m
%
% Function to plot all kernels, variance explained scatterplots, and
%  measured vs. predicted responses for all kernels when given kernels
%  struct and the corresponding kernelParams and kernelResampInputs structs
%  generated by extractAllKernels().
% Calls plot1KernelPred()
% NOTE: for a fly where both the left and right cells are included, this
%  generates 24 separate figures
% 
% INPUTS:
%   kernels - kernels struct from extractAllKernels()
%   kernelParams - kernels parameter struct from extractAllKernels()
%   kernelResampInputs - kernels inputs, resampled at kernel time scale,
%       from extractAllKernels()
% OUTPUTS:
%   none, but produces figures
%
% CREATED: 5/24/19
% UPDATED: 5/24/19
%

function plotAllKernelPred(kernels, kernelParams, kernelResampInputs)
    
    % get all fields in kernels struct (which kernels computed)
    kernelFields = fieldnames(kernels);
    
    % loop through all structs in kernels array (left, right, sum, diff)
    for i = 1:length(kernels)
        % which data
        switch i
            case 1
                dataName = 'Left';
            case 2
                dataName = 'Right';
            case 3
                dataName = 'Sum';
            case 4
                dataName = 'Diff';
        end
        
        % loop through all fields in kernels (for which kernels computed)
        for j = 1:length(kernelFields)
            
            % whether forward or reverse kernel, by f or r as first letter;
            %  use to determine which of inputs is actual response, whether
            %  to reverse kernel time axis
            switch kernelFields{j}(1)
                case 'f'
                    actResp = kernelResampInputs.dFF(i,:);
                    inv = 0; % no time axis reversal
                    
                case 'r'
                    if (strcmpi(kernelFields{j}, 'rFwdVel'))
                        actResp = kernelResampInputs.fwdVel;
                    elseif (strcmpi(kernelFields{j}, 'rYawVel'))
                        actResp = kernelResampInputs.yawAngVel;
                    elseif (strcmpi(kernelFields{j}, 'rYawSpd'))
                        actResp = kernelResampInputs.yawAngSpd;
                    end
                    inv = 1; % reverse kernel time axis
            end
            
            % plot title
            
            ttl = sprintf('%s %s', dataName, kernelFields{j});
            
            % produce plot
            plot1KernelPred(kernels(i).(kernelFields{j}).kernel, ...
                kernelParams.t, inv, ...
                kernels(i).(kernelFields{j}).predResp, actResp, ...
                    kernelResampInputs.t, ...
                    kernels(i).(kernelFields{j}).fitObj, ...
                    kernels(i).(kernelFields{j}).varExpl, ttl);
        end
    end

end

% plot1KernelPred.m
%
% Function that plots the kernel, the scatterplot of predicted vs. measured
%  responses with best fit line, and predicted and measured responses for 1
%  kernel.
% NOTE: kernel and responses should be on same time basis
%
% INPUTS:
%   kernel 
%   lags - time points for kernel
%   inv - flag for whether to invert kernel relative to lags
%   predResp - predicted response
%   measResp - measured response
%   respT - time points for predicted and measured responses
%   fitobj - fit object for best fit line
%   varExpl - variance explained
%   ttl - plot title
%
% OUTPUTS:
%   f - figure handle
%   also, generates plot
%
% CREATED: 5/23/19
% UPDATED: 5/23/19
%

function f = plot1KernelPred(kernel, lags, inv, predResp, measResp, ...
    respT, fitObj, varExpl, ttl)
    
    if inv
        kernel = fliplr(kernel);
    end

    f = figure;
    
    % plot kernel
    subplot(2,2,1);
    plot(lags, kernel);
    xlabel('Time (s)');
    
    % plot scatterplot of predicted vs. measured and best fit line
    subplot(2,2,2);
    scatter(predResp, measResp, 'Marker', '.', 'MarkerFaceAlpha',0.05, ...
        'MarkerEdgeAlpha', 0.05);
    hold on;
    plot(fitObj);
    legend('hide');
    xlabel('Predicted');
    ylabel('Measured');
    title(sprintf('Variance Explained = %.2f', varExpl));
    
    % plot actual predicted and measured responses
    subplot(2,2,[3,4]);
    plot(respT, measResp, 'b'); % measured response in blue
    hold on;
    plot(respT, predResp, 'r'); % predicted response in red
    xlabel('Time (s)');
    legend('Measured', 'Predicted');
    
    % title for whole plot
    suptitle(ttl);
end
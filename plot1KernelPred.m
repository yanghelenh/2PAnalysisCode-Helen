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
%   predResp - predicted response
%   measResp - measured response
%   respT - time points for predicted and measured responses
%   fitobj - fit object for best fit line
%   varExpl - variance explained
%   title - plot title
%
% OUTPUTS:
%   f - figure handle
%   also, generates plot
%
% CREATED: 5/23/19
% UPDATED: 5/23/19
%

function f = plot1KernelPred(kernel, lags, predResp, measResp, respT, ...
    fitobj, varExpl, title)
    
    f = figure;
    
    % plot kernel
    subplot(2,2,1);
    plot(lags, kernel);
    xlabel('Time (s)');
    
    % plot scatterplot of predicted vs. measured and best fit line
    subplot(2,2,2);
    plot(fitobj, predResp, measResp);
    legend('hide');
    xlabel('Predicted');
    ylabel('Measured');
    title(sprintf('Variance Explained = %.2f', varExpl));
    
    % plot actual predicted and measured responses
    subplot(2,2,[3,4]);
    hMeas = plot(respT, measResp, 'b'); % measured response in blue
    hold on;
    hPred = plot(respT, predResp, 'r'); % predicted response in red
    xlabel('Time (s)');
    ylabel('dF/F');
    legend({hMeas, hPred}, 'Measured', 'Predicted');
    
    % title for whole plot
    suptitle(title);
end
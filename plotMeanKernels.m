% plotMeanKernels.m
%
% Function to plot all kernels, where those relating the same variables in
%  the same direction are plotted on the same subplot. I.e. the yaw
%  velocity forward kernel for left, right, sum, and diff. Also shows
%  number of flies. Produces 1 plot, with subplots.
% Requires the structs generated by extractKernels()
%
% INPUT:
%   kernels - kernels struct from extractKernels()
%   kernelParams - kernel parameters struct
%   allYLims - cell array of yLims for all subplots, if [], let matlab
%       scale
%   ttl - title for whole plot
%   withSEM - boolean for whether to include error shading
%   withIndiv - boolean for whether to include individual flies
%
% OUTPUT:
%   f - handle to figure
%   also produces plot
%
% CREATED: 8/28/19 - HHY
% UPDATED: 8/28/19 - HHY
%

function f = plotMeanKernels(kernels, kernelParams, allYLims, ttl, ...
    withSEM, withIndiv)
    % get all fields in kernels struct (which kernels computed)
    kFN = fieldnames(kernels);

    numSubplots = length(kFN);
    subplotRows = 2;
    subplotCols = ceil(numSubplots/subplotRows);
    
    % check if user provided yLims
    useYLims = ~isempty(allYLims);
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    cmap = colormap('lines'); % get colormap
    
    for i = 1:numSubplots
        whichSubplot = i;
        
        % flip axis for reverse kernels
        if (kFN{i}(1) == 'r')
            lags = -1 * kernelParams.t;
        else
            lags = kernelParams.t;
        end
            
        
        subplot(subplotRows, subplotCols, whichSubplot);
        hold on;
        
        cFN = fieldnames(kernels.(kFN{i}));
        
        
        for j = 1:length(cFN)
            if withSEM
                lineHand(j) = plot_err_patch_v2(lags, ...
                    kernels.(kFN{i}).(cFN{j}).meanKernel, ...
                    kernels.(kFN{i}).(cFN{j}).sem, cmap(j,:) * 0.8, ...
                    cmap(j,:));
            else
                lineHand(j) = plot(lags, ...
                    kernels.(kFN{i}).(cFN{j}).meanKernel, ...
                    'Color', cmap(j,:), 'LineWidth', 1.5);
            end
            
            hold on;
            
            if withIndiv
                plot(lags, ...
                    kernels.(kFN{i}).(cFN{j}).allKernels, ...
                    'Color', cmap(j,:), 'LineWidth', 0.3);
            end
            
            legendTxt{j} = sprintf('%s %.2f (%d)', cFN{j}, ...
                mean(kernels.(kFN{i}).(cFN{j}).varExpl), ...
                kernels.(kFN{i}).(cFN{j}).numFlies);
                 
        end
        
        xlabel('Time (s)');
        xlim([-1*kernelParams.winLen kernelParams.winLen]);

        if (useYLims)
            y = allYLims{i};
        else
            y = ylim;
        end
        
        % x-axis line
        line([-1*kernelParams.winLen kernelParams.winLen], [0,0],...
            'Color', 'black');
        
        % y-axis line
        line([0,0], y, 'Color', 'black');
        ylim(y);
        
        title(kFN{i});
        legend(lineHand, legendTxt);
       
        
    end
    
    suptitle(ttl);
    
end
% plotLRrunCorr.m
%
% Function for looking at correlation in activity between left and right
%  cells in calcium imaging.
% Prompts user to select pData file through GUI and outputs graph showing
%  two dF/F traces, FicTrac yaw and forward velocities, and running
%  correlation. Also shows overall correlation b/w two cells.
%
% INPUTS:
%   corrWin - length of window, in seconds, over which to compute running
%       correlation
%   prompts for pData file through GUI
%
% OUTPUTS:
%   none, but produces graph
%
% CREATED: 9/10/21 - HHY
%
% UPDATED:
%   9/10/21 - HHY
%
function plotLRrunCorr(corrWin)

    % prompt user to select pData file
    disp('Select pData file');
    [pDatFile, pDatPath] = uigetfile('*.mat', 'Select a pData file', ...
        pDataPath());
    
    % load data from pData file
    load(fullfile(pDatPath, pDatFile), 'img', 'fictrac', 'name');
    
    
    % compute running correlation b/w left and right cells
    
    % get size of window, as odd number of frames
    ifi = median(diff(img.t)); % interframe interval
    corrWinFr = floor(corrWin / ifi);
    % if even, add 1
    if ~(mod(corrWinFr,2))
        corrWinFr = corrWinFr + 1;
    end
    halfWinSize = floor(corrWinFr / 2);
    
    
    % preallocate
    runCorr = zeros(size(img.dFF.left));
    
    % loop through all frames; window centered on that imaging frame
    for i = 1:length(img.dFF.left)
        startInd = i - halfWinSize; % get start index of window
        % window start would be before start of trial, just set to first
        %  frame
        if (startInd < 1)
            startInd = 1;
        end
        
        endInd = i + halfWinSize; % get end index of window
        % window start would be after end of trial, set to end
        if (endInd > length(img.dFF.left))
            endInd = length(img.dFF.left);
        end
        
        % dF/F snippet for left cell
        leftSeg = img.dFF.left(startInd:endInd);
        % dF/F snippet for right cell
        rightSeg = img.dFF.right(startInd:endInd);
        
        % get correlation b/w these two snippets
        thisCorrCoef = corrcoef(leftSeg,rightSeg);
        runCorr(i) = thisCorrCoef(1,2);   
    end
    
    % overall correlation
    overallCorr = corrcoef(img.dFF.left,img.dFF.right);
    overallCorr = overallCorr(1,2);
    
    
    % plot
    
    figure;
    
    % dF/F for left and right cells
    ax1 = subplot(4,1,1);
    plot(img.t, img.filtDFF.left);
    hold on;
    plot(img.t, img.filtDFF.right);
    yScaleMin = prctile([img.filtDFF.left img.filtDFF.right], 5) * 1.2;
    yScaleMax = prctile([img.filtDFF.left img.filtDFF.right], 95) * 1.2;
    ylim([yScaleMin yScaleMax]);
    legend({'left','right'});
    title('dF/F');
    
    % plot running correlation
    ax2 = subplot(4,1,2);
    plot(img.t, runCorr);
    hold on;
    line([img.t(1) img.t(end)], [0 0], 'Color','black');
    title(sprintf('Running Correlation, overall = %0.3f',overallCorr));
    
    % plot FicTrac yaw
    ax3 = subplot(4,1,3);
    plot(fictrac.t, fictrac.yawAngVel);
    hold on;
    line([fictrac.t(1) fictrac.t(end)], [0 0], 'Color','black');
    ylim([-500 500]);
    title('FicTrac Yaw Velocity');
    
    % plot FicTrac fwd
    ax4 = subplot(4,1,4);
    plot(fictrac.t, fictrac.fwdVel);
    hold on;
    line([fictrac.t(1) fictrac.t(end)], [0 0], 'Color','black');
    ylim([-5 20]);
    title('FicTrac Forward Velocity');
    
    % link axes
    linkaxes([ax1 ax2 ax3 ax4], 'x');
    
end
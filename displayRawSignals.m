% displayRawSignals.m
%
% Function that plots dF/F over time for each ROI, as well as raw
%  fluorescence before background subtraction as separate plot. Displays
%  average image with ROIs overlaid
%
% INPUT:
%   frameTimes - time at which each frame occured
%   avSignal - raw fluorescence before background subtraction
%   bksSignal - background subtracted signal
%   roiMasks - masks for all ROIs
%   avImg - mean aligned image
%
% OUTPUT:
%   none but generates figures as side effects
%
% CREATED: 12/10/18, but adapted from grad school version
% UPDATED: 12/10/18
%

function displayRawSignals(frameTimes, avSignal, bksSignal, roiMasks,...
    avImg)

    nMasks = length(roiMasks);

    % plot signals 
    cm = colormap('lines');

    % plots intensity over time, no background subtraction (avSignal)
    for i = 1:nMasks
        if (i>64)
            cmInd = mod(i,65)+1;
        else
            cmInd = i;
        end
        plot(frameTimes, avSignal(i,:), 'color', cm(cmInd,:)); 
        hold on;
    end
    xlabel('Time (sec)');
    title('Intensity in ROIs - before background substraction');

    % plot dF/F
    figure; 
    for i = 1:nMasks
        if (i>64)
            cmInd = mod(i,65)+1;
        else
            cmInd = i;
        end
        % ***NOTE: this dF/F calculation is for visualization purposes. It
        %  assumes the baseline F is the mean intensity of the bkgnd-
        %  subtracted average image. Depending on the stimulus, the 
        %  baseline F may not be the mean! 
        % dF/F for each ROI
        trace = bksSignal(i,:)/mean(bksSignal(i,:)) + 1*(i-1); 
        plot(frameTimes, trace, 'color', cm(cmInd,:), 'linewidth', 2);
        hold on;
    end

    % set axis limits (x min, x max, y min, y max)
    axis([0 frameTimes(end) 0 i+1]);

    xlabel('time (sec)'); 
    ylabel('ROI #');
    title('dF/F');

    % Create a colored map of ROIs
    xPixels = size(roiMasks{1}, 1);
    yPixels = size(roiMasks{1}, 2);
    CMask = zeros(xPixels, yPixels, 3);
    
    for i = 1:nMasks
        if (i>64)
            cmInd = mod(i,65)+1;
        else
            cmInd = i;
        end
        curColor = cm(cmInd,:);
        curMask = cat(3, curColor(1).*roiMasks{i}, ...
            curColor(2).*roiMasks{i}, curColor(3).*roiMasks{i});
        CMask = CMask + curMask;
    end

    figure;
    imshow(avImg,[]); 
    hold on;
    h = imshow(CMask);
    set(h,'AlphaData',0.5);
    title('ROI masks');

end
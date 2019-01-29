% plotRawImgFicTracData.m
%
% Function to plot raw traces of imaging data (dF/F) for all ROIs and 
%  FicTrac data (yaw angular velocity and forward velocity), each in its 
%  own separate subplot
%
% INPUTS:
%   frameTimes - frame times of imaging data
%   bksSignal - background subtracted fluorescence signal for all ROIs,
%       from uSelectROIs
%   fictracTimes - timing for FicTrac data
%   fwdVel - forward velocity, in mm/sec
%   yawAngVel - rotational angular velocity, in degrees/sec
%
% OUTPUTS:
%   none, but generates plot as side effect
%
% CREATED: 12/12/18 HHY
% UPDATED: 12/12/18 HHY
%

function plotRawImgFicTracData(frameTimes, bksSignal, fictracTimes, ...
    fwdVel, yawAngVel, tRange)
    
    numROIs = size(bksSignal, 1); 
    
    numSubplots = numROIs + 2; % number of ROIs plus 2 for FicTrac data
    
    dFFs = zeros(size(bksSignal));
    % compute dF/F for all ROIs
    for i = 1:numROIs
        dFFTemp = computeDFF(bksSignal(i,:), frameTimes, ...
            bksSignal(i,:), frameTimes);
        dFFs(i,:) = dFFTemp';
    end
    
    % determine ylim on dF/F subplots
    yMax = max(max(dFFs));
    yMin = min(min(dFFs));
    
    % determine xlim on all plots
    xMax = max([frameTimes(end), fictracTimes(end)]);
%     xMin = min([frameTimes(1), fictracTimes(1)]);

    % determine starting xLim on all plots
    tLims = [0 tRange];
    
    f = figure;
    
    % plot ROIs dF/F
    for i = 1:numROIs
        subplotHandles{i} = subplot(numSubplots, 1, i);
        plot(frameTimes, dFFs(i,:));
%         ylim([yMin, yMax]);
        ylim([yMin, 1]);
        xlim(tLims);
%         xlim([xMin, xMax]);
        xlabel('Time (sec)');
        ylabel('dF/F');
        title(sprintf('ROI %d', i)); 
    end
    
    % plot FicTrac data
    subplotHandles{numROIs+1} = subplot(numSubplots, 1, numROIs + 1);
    plot(fictracTimes, fwdVel);
%     xlim([xMin, xMax]);
    xlim(tLims);
    ylim([-10 30]);
    xlabel('Time (sec)');
    ylabel('mm/sec');
    title('Forward Velocity');
    
    subplotHandles{numROIs+2} = subplot(numSubplots, 1, numROIs + 2);
    plot(fictracTimes, yawAngVel);
%     xlim([xMin, xMax]);
    ylim([-500 500]);
    xlim(tLims);
    xlabel('Time (sec)');
    ylabel('deg/sec');
    title('Rotational Velocity');
    
    tSlider = uicontrol(f, 'Style', 'slider', 'Position', [20 10 400 20]);
    tSlider.Value = 0;
    tSlider.Callback = @updateTLim;
    
    function updateTLim(src, event)
        for j = 1:length(subplotHandles)
            xlim(subplotHandles{j}, ...
                [tSlider.Value * xMax, tSlider.Value * xMax + tRange]);
        end
    end
    

end
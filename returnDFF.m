% returnDFF.m
%
% Function that converts background-subtracted fluorescence signal into
%  dF/F (for 1 channel data) or dR/R (for 2 channel data) for a single
%  trial folder. When there is only 1 ROI, prompts user for whether it's
%  left or right cell. When there are 2 ROIs, also returns sum and
%  difference of dF/F (difference as right - left).
% Returns dF/F, timing info, number of channels in data, the average image,
%  and the roiMasks.
% Note: assumes that every imaging frame should be used for the baseline
%  calculation. If this is not the case, update this function or write a
%  new one.
% Adaptation of computeDFF_trial.m, only instead of saving the data to a
%  pData file, returns it.
%
% INPUTS:
%   trialPath - full path to trial folder to compute dF/F
%
% OUTPUTS:
%   dFF - struct of dF/F (fields: left, right, sum, diff or roi# for when
%       number of ROIs > 2)
%   t - frame timing (as average b/w frame start and frame end times
%   numChannels - number of channels in imaging
%   numROIs - number of ROIs
%   avgImg - average image (for both channels if 2 channel data), passed
%       directly from imDat.mat
%   roiMasks - masks for ROIs, passed directly from imDat.mat
%
% CREATED: 8/27/19 - HHY
% 
% UPDATED: 8/27/19 - HHY
%

function [dFF, t, numChannels, numROIs, avgImg, roiMasks] = returnDFF(...
    trialPath)

    curDir = pwd;
    
    cd(trialPath)

    % load data
    load('imDat.mat', 'bksSignal', 'frameStartTimes', 'frameEndTimes',...
        'roiMasks', 'meanImgAligned');
    
    % rename avgImg
    avgImg = meanImgAligned;
    
    % compute t as mean b/w frame start and end times
    t = mean([frameStartTimes; frameEndTimes], 1);

    % compute dF/F or dR/R

    % determine number of channels 
    bksSignalFields = fieldnames(bksSignal);
    numChannels = length(bksSignalFields);
    
    % determine number of ROIs
    numROIs = size(bksSignal.ch1, 1);

    % single-channel data
    if (numChannels == 1)
        
        % pre-allocate
        dFFs = zeros(size(bksSignal.ch1));

        % compute dF/F
        for i = 1:numROIs
            ch1Traces = bksSignal.ch1;
            dFFs(i,:) = computeDFF(ch1Traces(i,:), t, ch1Traces(i,:), t)';
        end

    % 2-channel data
    elseif (numChannels == 2)

        % number of ROIs
        numROIs = size(bksSignal.ch1, 1); 
        
        % pre-allocate
        dFFs = zeros(size(bksSignal.ch1));

        % compute dR/R (though named dF/F)
        ch1Traces = bksSignal.ch1;
        ch2Traces = bksSignal.ch2;
        for i = 1:numROIs
            ch1dFFTemp = computeDFF(ch1Traces(i,:), ...
                t, ch1Traces(i,:), t) + 1;
            ch2dFFTemp = computeDFF(ch2Traces(i,:), ...
                t, ch2Traces(i,:), t) + 1;
            ch1DivCh2Temp = ch1dFFTemp ./ ch2dFFTemp - 1; 
            dFFs(i,:) = ch1DivCh2Temp';
        end
    % strange output, not 1 or 2 channel data
    else
        fprintf('%s \n does not have bksSignal with 1 or 2 channels \n', ...
            trialPath);
        dFF = []; % return empty dF/F
        cd(curDir)
        return;
    end
    
    % convert dFFs to dFF struct, label left and right cell appropriately
    if (numROIs == 1)
        disp(trialPath);
        prompt = 'Is the single ROI the left or right cell? (L/R) ';
        ui = input(prompt,'s');
        if (strcmpi(ui, 'L'))
            dFF.left = dFFs;
        elseif (strcmpi(ui, 'R'))
            dFF.right = dFFs;
        else
            dFF = [];
        end
    elseif (numROIs == 2)
        dFF.left = dFFs(1,:);
        dFF.right = dFFs(2,:);
        dFF.sum = dFFs(1,:) + dFFs(2,:);
        dFF.diff = dFFs(2,:) - dFFs(1,:);
    else % if there are more than 2 ROIs, just name them as struct fields
        for i = 1:numROIs
            fieldName = sprintf('roi%d', i);
            dFF.(fieldName) = dFFs(i,:);
        end
    end

    cd(curDir);
end
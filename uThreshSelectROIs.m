% uThreshSelectROIs.m
%
% Function for the user to select the ROIs and background threshold for a 
%  single time series. Calculates raw signal (mean pixel intensity) in each
%  ROI and performs background subtraction.
% ROI selection by pixels with mean pixel intensity relative to rest of
%  image above threshold and within user-drawn outline.
%
% INPUTS:
%   none, but prompts user for trial folder with preprocessed imaging data
%
% OUTPUTS:
%   none, but saves data for ROIs back into same imDat.mat file
%
% CREATED: 2/13/19 HHY
% UPDATED: 3/7/19 HHY
%

function uThreshSelectROIs()

    close all % close all open figures
    clc
    
    % prompt for trial folder
    disp('Select a trial folder to preprocess.');
    uTrialPath = uigetdir;
    
    fprintf('Processing %s \n', uTrialPath);
    
    % check that selected folder has imDat.mat file
    trialPathFiles = dir(uTrialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasImDat = sum(strcmp(trialPathFileNames, 'imDat.mat'));
    
    if (hasImDat) % has appropriate file
        curDir = pwd;
        cd(uTrialPath);
        
        % get variables in imDat
        imDatStrct = whos('-file', 'imDat.mat');
        imDatVars = extractfield(imDatStrct, 'name');
        
        % check if background subtracted signal exists (i.e. ROIs selected
        %  before)
        roisSelected = sum(strcmp(imDatVars, 'bksSignal'));
        
        % if background subtracted signal exists, prompt user if they want
        %  to overwrite
        if (roisSelected) 
            prompt = ['ROIs have already been selected for this trial.' ... 
                ' Overwrite? (Y/N) '];
            ui = input(prompt,'s');
            if (~strcmpi(ui, 'Y'))
                % stop running this function. don't overwrite 
                disp('Ending uThreshSelectROIs. Nothing overwritten');
                cd(curDir); % return to previous directory
                return;
            end
        end
        
        % load imDat.mat contents, only needed variables
        load('imDat.mat', 'alignedSeries', 'frameStartTimes', ...
            'meanImgAligned');

        % combination draw and thresh to get masks for ROIs - always on ch1
        roiMasks = threshDrawROIs(meanImgAligned.ch1);

        % user determines background (by setting threshold on dF/F
        % of image across space
        bkMask = selectBackground(meanImgAligned.ch1);

        numChannels = numel(fieldnames(meanImgAligned));
        
        % get raw signal for each ROI for each channel and plot
        for i = 1:numChannels
            
            channelName = ['ch' num2str(i)];
            
            % get average signal and background subtracted signal for
            %  each ROI
            [avSignal.(channelName), bksSignal.(channelName)] = ...
                getRawSignals(alignedSeries.(channelName), ...
                roiMasks, bkMask);
            
            % plot dF/F for each ROI
            displayRawSignals(frameStartTimes, avSignal.(channelName), ...
                bksSignal.(channelName), roiMasks, ...
                meanImgAligned.(channelName))
        end

        % save new variables into imDat.mat
        save('imDat.mat', 'roiMasks', 'bkMask', 'avSignal', ...
            'bksSignal', '-append');
        disp('Saved!');
        
        cd(curDir); % return to previous directory
        
    else
        disp('Selected trial folder does not contain imDat.mat file');
        return;
    end
        
end
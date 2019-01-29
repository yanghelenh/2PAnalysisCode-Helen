% uSelectROIs.m
%
% Function for the user to select the ROIs and background threshold for a 
%  single time series. Calculates raw signal (mean pixel intensity) in each
%  ROI and performs background subtraction.
%
% INPUTS:
%   none, but prompts user for trial folder with preprocessed imaging data
%
% OUTPUTS:
%   none, but saves data for ROIs back into same imDat.mat file
%
% CREATED: 12/10/18 HHY
% UPDATED: 12/10/18 HHY
%

function uSelectROIs()

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
        % load imDat.mat contents
        load('imDat.mat');
        
        % if background subtracted signal exists (i.e. uSelectROIs already
        %  run on this trial folder before)
        if (exist('bksSignal', 'var')) 
            prompt = ['ROIs have already been selected for this trial.' ... 
                ' Overwrite? (Y/N) '];
            ui = input(prompt,'s');
            if (~strcmpi(ui, 'Y'))
                % stop running this function. don't overwrite 
                disp('Ending uSelectROIs. Nothing overwritten');
                cd(curDir); % return to previous directory
                return;
            end
        end
        
        
        % user draws ROIs, get masks for ROIs
        roiMasks = drawROIs(meanImageAligned);

        % user determines background (by setting threshold on dF/F
        % of image across space
        bkMask = selectBackground(meanImageAligned);

        % get average signal and background subtracted signal for
        %  each ROI
        [avSignal, bksSignal] = getRawSignals(alignedSeries, ...
            roiMasks, bkMask);

        % plot dF/F for each ROI
        displayRawSignals(frameStartTimes, avSignal, bksSignal, ...
            roiMasks, meanImageAligned)

        % save new variables into imDat.mat
        save('imDat.mat', 'unalignedSeries', 'alignedSeries', ...
            'meanImageAligned', 'trialPath', 'tifFileName', ...
            'frameStartTimes', 'frameEndTimes', 'frameTimingError', ...
            'roiMasks', 'bkMask', 'avSignal', 'bksSignal', '-v7.3');
        disp('Saved!');
        
        cd(curDir); % return to previous directory
        
    else
        disp('Selected trial folder does not contain imDat.mat file');
        return;
    end
        
end
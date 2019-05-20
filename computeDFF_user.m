% computeDFF_user.m
%
% Function that converts background-subtracted fluorescence signal into
%  dF/F (for 1 channel data) or dR/R (for 2 channel data) for user selected
%  trial folder. Saves processed data in pData file in same folder.
% Note: assumes that every imaging frame should be used for the baseline
%  calculation. If this is not the case, update this function or write a
%  new one.
%
% INPUTS:
%   None, but prompts user to select trial folder to compute dF/F on.
%
% OUTPUTS:
%   None, but creates pData.mat file to save dF/F, timing, and # of
%     channels info
%
% CREATED: 3/7/19 HHY
% UPDATED: 3/25/19 HHY
%

function computeDFF_user()

    % ask user to select trial folder
    disp('Select a trial folder analyze.');
    uTrialPath = uigetdir;
    curDir = pwd;
    cd(uTrialPath)
    
    fprintf('Computing dF/F on %s \n', uTrialPath);
    
    % check that selected folder has imDat.mat file
    trialPathFiles = dir(uTrialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasImDat = sum(strcmp(trialPathFileNames, 'imDat.mat'));
    
    if (hasImDat)
        % get variables in imDat
        imDatStrct = whos('-file', 'imDat.mat');
        imDatVars = extractfield(imDatStrct, 'name');
        
        % check if background subtracted signal exists (i.e. ROIs selected
        %  before)
        roisSelected = sum(strcmp(imDatVars, 'bksSignal'));
        if (roisSelected)
            
            % load data
            load('imDat.mat', 'bksSignal', 'frameStartTimes');
            
            % compute dF/F or dR/R
            
            % determine number of channels 
            bksSignalFields = fieldnames(bksSignal);
            numChannels = length(bksSignalFields);

            % single-channel data
            if (numChannels == 1)
                
                % number of ROIs
                numROIs = size(bksSignal.ch1, 1);
                
                % pre-allocate
                dFFs = zeros(size(bksSignal.ch1));
                
                % compute dF/F
                for i = 1:numROIs
                    ch1Traces = bksSignal.ch1;
                    dFFs(i,:) = computeDFF(ch1Traces(i,:), ...
                        frameStartTimes, ch1Traces(i,:), ...
                        frameStartTimes)';
                end
            
            % 2-channel data
            elseif (numChannels == 2)
            
                % number of ROIs
                numROIs = size(bksSignal.ch1, 1);               
                 
                % compute dR/R (though named dF/F)
                ch1Traces = bksSignal.ch1;
                ch2Traces = bksSignal.ch2;
                for i = 1:numROIs
                    ch1dFFTemp = computeDFF(ch1Traces(i,:), ...
                        frameStartTimes, ch1Traces(i,:), ...
                        frameStartTimes) + 1;
                    ch2dFFTemp = computeDFF(ch2Traces(i,:), ...
                        frameStartTimes, ch2Traces(i,:), ...
                        frameStartTimes) + 1;
                    ch1DivCh2Temp = ch1dFFTemp ./ ch2dFFTemp - 1; 
                    dFFs(i,:) = ch1DivCh2Temp';
                end
            
            else
                disp('bksSignal does not have 1 or 2 channels.');
                cd(curDir)
                return;
            end
            
            % save data
            
            % check if pData.mat file already exists
            hasPData = sum(strcmp(trialPathFileNames, 'pData.mat'));
            
            % don't overwrite whole pData file, just the dF/F variables if
            %  pData.mat already exists
            if (hasPData)
                save('pData.mat', 'dFFs', 'numChannels', ...
                    'frameStartTimes', '-append');
                disp('Saved!');
            else
                save('pData.mat', 'dFFs', 'numChannels', ...
                    'frameStartTimes', '-v7.3');
                disp('Saved!');
            end
            
            cd(curDir);

        else
            disp('ROIs have not been selected for this trial');
            cd(curDir);
            return;
        end
        
    else
        disp('Selected trial folder does not contain imDat.mat file');
        cd(curDir);
        return;
    end
end
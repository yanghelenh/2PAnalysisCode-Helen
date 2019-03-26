% analyzeBasicFWA_fictrac.m
%
% Function to compute fluorescence weighted average for each ROI on the
%  FicTrac parameters. Saves the output in a pData file in the same folder.
%  
% Note that this is the basic version of this function, such that it does
%  not account for temporal correlations in the FicTrac outputs. It also
%  doesn't account for temporal correlations in the ROI responses,
%  especially, for example, those imposed by the indicator. 
%
% Currently (3/7/19), the code isn't written to exclude times when FicTrac
%  dropped and therefore the FicTrac output values are meaningless. Update
%  this.
%
% INPUTS:
%   timeWindow - length of extracted FWA on either side of 0 (i.e. total
%       'filter' length is 2x timewindow)
%   smoAvgWindow - length, in seconds, of moving average window used to
%       smooth FicTrac data
%   
%
% OUTPUTS:
%   none, but saves output in pData file in same folder and produces plots
%    of FWAs
%
% CREATED: 3/7/19 HHY
% UPDATED: 3/7/19 HHY
%

function analyzeBasicFWA_fictrac(timeWindow, smoAvgWindow)

    % ask user to select trial folder
    disp('Select a trial folder to analyze.');
    uTrialPath = uigetdir;
    curDir = pwd;
    cd(uTrialPath)
    
    fprintf('Analyzing %s \n', uTrialPath);
    
    % check that selected folder has imDat.mat and fictracDat.mat files
    trialPathFiles = dir(uTrialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasImDat = sum(strcmp(trialPathFileNames, 'imDat.mat'));
    hasFictracDat = sum(strcmp(trialPathFileNames, 'fictracDat.mat'));
    
    % if both files exist in this trial folder
    if (hasImDat && hasFictracDat)
        % check if imDat.mat has had ROIs selected
        % get variables in imDat
        imDatStrct = whos('-file', 'imDat.mat');
        imDatVars = extractfield(imDatStrct, 'name');
        
        % check if background subtracted signal exists (i.e. ROIs selected
        %  before)
        roisSelected = sum(strcmp(imDatVars, 'bksSignal'));
        
        % only process if ROIs have been selected
        if (roisSelected)
            % load data
            load('fictracDat.mat', 'fwdVel', 'yawAngVel', 't');
            fictracTimes = t;
            load('imDat.mat', 'bksSignal', 'frameStartTimes');
            frameTimes = frameStartTimes;
            
            numROIs = size(bksSignal.ch1, 1); 
                
            % compute dF/F
            
            
        else
            disp('ROIs have not been selected for this trial');
            cd(curDir);
            return;
        end
        
    % save pData file
    % if pData file already exists, use -append
    % otherwise, save new file
    
    else
        disp('Selected trial folder does not contain necessary files');
        cd(curDir);
        return;
    end
end
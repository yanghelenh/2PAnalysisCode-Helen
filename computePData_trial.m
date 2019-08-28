% computePData_trial.m
%
% Function to compute all processed data (currently, for imaging and
%  FicTrac data) for a single trial and save it into a pData.mat file in
%  the same trial folder.
% Preprocessing, ROI selection, and FicTrac dropped frames selection must
%  have been completed
%
% INPUT:
%   trialPath - full path of a single trial folder
%
% OUTPUT:
%   none, but generates pData.mat file in trial folder containing:
%   
%   img - struct containing all img data
%       dFF - dF/F for all ROIs as struct; if numROIs == 2, includes sum 
%           b/w 2 ROIs and difference (ROI2 - ROI1)
%       numChannels - number of imaging channels (1 or 2)
%       t - frame times (average b/w frame start and end times)
%       filtDFF - gaussian process smoothed dFF
%       numROIs - number of ROIs
%       avgImg - average image across whole time series
%       roiMasks - masks for ROI
%       filtParams - struct of parameters for smoothing dFF
%           padLen - padding length for convolving with Gaussian kernel
%           sigma - standard deviation for Gaussian kernel
%           gpRatio - ratio for gaussian process smoothing
%           priorStdDev - prior standard deviationm, gaussian process
%               smoothing
%   fictrac - struct containing downsampled and smoothed FicTrac data
%       dsf - downsampling factor
%       fwdCumPos
%       fwdVel
%       slideCumPos
%       slideVel
%       yawAngCumPos
%       yawAngPosWrap
%       yawAngVel
%       xPos
%       yPos
%       yawAngSpd
%       totAngSpd
%       totSpd
%       t
%       dropInd
%       filtParams - struct of parameters for smoothing all fictrac data
%           padLen - padding length for convolving with Gaussian kernel
%           sigmaPos - standard deviation for Gaussian kernel for position
%           sigmaVel - standard deviation for Gaussian kernel for velocity
%       degPerMM - degrees per millimeter, conversion factor
%       mmPerDeg - millimeters per degree, conversion factor
%   exptCond - experiment type, as string
%   name - struct containing name info for trial, all strings
%       dateName
%       flyName
%       fovName
%       trialName
%       exptName - full experiment name DATE_fly##_fov##_trial##
%
% CREATED: 8/27/19 - HHY
% UPDATED: 8/27/19 - HHY

function computePData_trial(trialPath)

    % CONSTANTS
    img.filtParams.padLen = int32(10);
    img.filtParams.sigma = int32(1);
    img.filtParams.gpRatio = 0.2;
    img.filtParams.priorStdDev = 100;

    fictrac.dsf = 20; % downsample to 500 Hz;
    fictrac.filtParams.padLen = int32(100);
    fictrac.filtParams.sigmaPos = int32(50); % 100 ms
    fictrac.filtParams.sigmaVel = int32(25); % 50 ms
    
    BALL_DIAM = 6.46;
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    fictrac.mmPerDeg = circum / 360; % mm per degree of ball
    fictrac.degPerMM = 360 / circum; % deg per mm ball

    % save current directory, return to it at end
    curDir = pwd;

    cd(trialPath);
    
    % check that the the current directory has an imDat.mat file as well as
    %  a fictracDat.mat file
    trialPathFiles = dir(trialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasImDat = sum(strcmp(trialPathFileNames, 'imDat.mat'));
    hasFTDat = sum(strcmp(trialPathFileNames, 'fictracDat.mat'));
    
    if (hasImDat && hasFTDat)
        % check that ROI selection and FicTrac dropped frames selection has
        %  been performed
        if (isempty(who('-file', 'fictracDat.mat', 'dropInd')))
            fprintf('%s \nlacks FicTrac dropInd\n', trialPath);
            return;
        end
        if (isempty(who('-file', 'imDat.mat', 'bksSignal')))
            fprintf('%s \nhas not had ROIs selected\n', trialPath);
            return;
        end
    
        % compute dF/F
        [dFF, t, numChannels, numROIs, avgImg, roiMasks] = ...
            returnDFF(trialPath);

        % gaussian process smooth dF/F
        dFFfieldNames = fieldnames(dFF);
        % loop through all fields
        for i = 1:length(dFFfieldNames)
            dat = dFF.(dFFfieldNames{i});

            % run light smoothing by convolving with gaussian kernel
            lightPySmo = py.proc_utils.safe_interp_conv_smooth_ball_data(...
                dat, img.filtParams.padLen, img.filtParams.sigma);
    %         lightSmo = cell2mat(cell(lightPySmo.tolist()));

            % run gaussian process smoothing
            gpPySmo = py.gp_smooth.gp_smooth(lightPySmo, ...
                img.filtParams.gpRatio, img.filtParams.priorStdDev);
            gpSmo = cell2mat(cell(gpPySmo.tolist()));

            filtDFF.(dFFfieldNames{i}) = gpSmo;
        end

        % save imaging data into struct
        img.dFF = dFF;
        img.filtDFF = filtDFF;
        img.t = t;
        img.numChannels = numChannels;
        img.numROIs = numROIs;
        img.avgImg = avgImg;
        img.roiMasks = roiMasks;
        
        % downsample and smooth FicTrac data, updates fictrac struct
        fictrac = dsFiltFictrac(trialPath, fictrac);
        
        % get experimental condition
        load('userDaqDat.mat', 'exptCond');
        
        % make name struct
        name = getExptName(trialPath);
        
        % save data, overrides any previous pData file
        save('pData.mat', 'img', 'fictrac', 'exptCond', 'name', '-v7.3');
        
    else % when imDat.mat and fictracDat.mat files missing
        fprintf(['%s \n does not contain imDat.mat file' ...
            ' and/or fictracDat.mat \n'], trialPath);
        cd(curDir);
        return;
    end
    
    
    % return to original directory
    cd(curDir);
end

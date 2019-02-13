% preprocessImaging.m
%
% Function for preprocessing just the Scanimage imaging data. Called by
%  preprocess(). Reads in ScanImage .tif files. Aligns each time series. 
%  Saves the aligned series, unaligned series, mean image, and metadata 
%  output by ScanImage. Handles both one and multichannel data. Alignment
%  performed on ch1 (as captured by ScanImage) using NoRMCorre. All other
%  channels shifted the same amount.
%
% NOTE: assumes that trial folder has only one f*.tif file, which is the
%  Scanimage .tif file. All NoRMCorre parameters at top.
%
% INPUTS:
%   tifFile - struct of info about the .tif file, output of dir()
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqTime - vector of times corresponding to each sample point of daqData
%
% OUTPUTS:
%   None, but saves data in imMetaDat.mat and imDat.mat for metadata and
%    aligned series, unaligned series, mean image imaging data,
%    respectively, in pwd.
%
% CREATED: 10/24/18 HHY
% UPDATED: 1/29/19 HHY
%

function preprocessImaging(tifFile, daqData, daqTime)
    % NoRMCorre Parameters here:
    BIN_WIDTH = 200; 
    MAX_SHIFT = 15;
    US_FAC = 50;
    INIT_BATCH = 200;

    tifFileName = tifFile.name;
    
    % print to command line
    fprintf('Reading %s \n', tifFileName);
    
    % read imaging data (ScanImage Tiff)
    imgData = ScanImageTiffReader(tifFileName);

    % get metadata
    metadata = imgData.metadata();
    % convert metadata into matlab struct from string
    % end of useful metadata
    endInd = strfind(metadata, '"RoiGroups"') - 5;
    % convert metadata to SI struct
    evalc(metadata(1:endInd));

    % get descriptions
    descriptions = imgData.descriptions();
    
    % convert descriptions into struct array
    % get current variables 
    currVars = who();
    
    % add all variables used in next part of function as well
    currVars = [currVars; {'currVars'; 'i'; 'j'; 'ans'; 'allVars'; ...
        'newVars';'tempStrct';'descriptionsStrct'}];
    
    for i = 1:length(descriptions)
        % run string for descriptions; loads variables into workspace
        evalc(descriptions{i});
        
        % identify which variables were just loaded
        allVars = who(); % all variables
        % variables loaded by descriptions
        newVars = allVars(not(cellfun(@(x) sum(strcmp(x, currVars)), ...
            allVars)));
        
        % construct strct from descriptions data
        for j = 1:length(newVars)
            tempStrct.(newVars{j}) = nestedEvalin(newVars{j});
        end
        descriptionsStrct(i) = tempStrct;
    end
 
    % save metadata and descriptions
    save('imMetaDat.mat', 'SI', 'descriptionsStrct', '-v7.3');
    
    % determine number of channels in data
    numChannels = length(SI.hChannels.channelSave);

    % get imaging data and permute so x is rows, y is columns
    imgDataSeries = permute(imgData.data(),[2,1,3]);
    
    % segment imaging data into separate 3D matricies for each channel
    %  save in struct 
    for i = 1:numChannels
        % name the channel using index from scanimage, prefixed with 'ch'
        channelName = ['ch' num2str(SI.hChannels.channelSave(i))];
        
        % indicies for this channel
        ind = i:numChannels:size(imgDataSeries,3);
        
        % this channel's time series
        tempSeries = imgDataSeries(:,:,ind);
        % convert from int16 to double
        tempSeries = double(tempSeries);
        % set minimum of series to 0
        tempSeries = tempSeries - min(tempSeries(:));
        
        % save series
        unalignedSeries.(channelName) = tempSeries;
    end
    
    % get number of frames in imaging data
    numFrames = size(unalignedSeries.ch1, 3);

    % align time series using NoRMCorre
    %  align on ch1, shift any additional channels to match
    % set parameters
    optionsRigid = NoRMCorreSetParms('d1', size(unalignedSeries.ch1,1),...
        'd2', size(unalignedSeries.ch1,2), 'bin_width', BIN_WIDTH, ...
        'max_shift', MAX_SHIFT, 'us_fac', US_FAC, 'init_batch', ...
        INIT_BATCH);
    % run rigid motion correction
    [alignedSeries.ch1, shifts, ~, optionsRigid] = normcorre(...
        unalignedSeries.ch1, optionsRigid);
    % compute motion correction metrics, provided by NoRMCorre
    [corrCoeffs.ch1, meanImgAligned.ch1, ~] = motion_metrics(...
        alignedSeries.ch1, 10);
    % align all other channels (if present) to ch1
    if (numChannels > 1)
        for i = 2:numChannels
            channelName = ['ch' num2str(i)];

            % shift channel
            alignedSeries.(channelName) = apply_shifts(...
                unalignedSeries.(channelName), shifts, optionsRigid);
            % compute motion correction metrics on other channels
            [corrCoeffs.(channelName), meanImgAligned.(channelName), ~]...
                = motion_metrics(alignedSeries.(channelName), 10);        
        end
    end
    
    % extract frame times
    frameStarts = find(diff(daqData.scanimageFrameClock) > 0.1);
    frameEnds = find(diff(daqData.scanimageFrameClock) < -0.1);
    
    frameTimingError = 0; % flag for timing error
    % should capture all frame starts and ends, so equal number
    if (length(frameStarts) ~= length(frameEnds))
        fprintf(['Unexpected ScanImage Frame Clock values for %s. \n'...
            'Number of frame starts and ends not equal.'], tifFileName);
        frameTimingError = 1;
    % after stop signal in middle of frame, takes time to stop, extra
    %  transition in frame clock; to be expected
    elseif (length(frameStarts) - 1 == numFrames)
        frameStarts = frameStarts(1:(end-1));
        frameEnds = frameEnds(1:(end-1));
    % if number of frames as indicated by frameStarts is not equal to the
    %  number of frames in the time series and also isn't the number of
    %  frames in the time series plus 1 (as above), then something is wrong
    %  in getting frame timing or in imaging file
    elseif (length(frameStarts) ~= numFrames)
        fprintf(['Number of frames in time series does not match '...
            'ScanImage Frame Clock values for %s.'], tifFileName);
        frameTimingError = 1;
    end
    
    % convert frame time indicies to actual time
    frameStartTimes = daqTime(frameStarts + 1);
    frameEndTimes = daqTime(frameEnds + 1);
    
    % trial path is present working directory
    trialPath = pwd;

    % save time series data
    save('imDat.mat', 'unalignedSeries', 'alignedSeries', ...
        'meanImgAligned', 'corrCoeffs', 'trialPath', 'tifFileName', ...
        'frameStartTimes', 'frameEndTimes', 'frameTimingError', '-v7.3');
    
    disp('Imaging data saved!');
                
end

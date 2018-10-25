% preprocessImaging.m
%
% Function for preprocessing just the Scanimage imaging data. Called by
%  preprocess(). Reads in ScanImage .tif files. Aligns each time series. 
%  Saves the aligned series, unaligned series, mean image, and metadata 
%  output by ScanImage. 
%
% NOTE: assumes that trial folder has only one .tif file, which is the
%  Scanimage .tif file. Also, note numRefFrames parameter for alignment at
%  top of function.
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
% UPDATED: 10/25/18 HHY
%

function preprocessImaging(tifFile, daqData, daqTime)

    numRefFrames = 30; % default to 30 reference frames for alignment

    tifFileName = tifFile.name;

    % read imaging data (ScanImage Tiff)
    imgData = ScanImageTiffReader(tifFileName);

    % get time series data
    fprintf('Reading %s \n', tifFileName);

    unalignedSeries = imgData.data();
    
    % number of frames in time series
    numFrames = size(unalignedSeries, 3);

    % get metadata
    metadata = imgData.metadata();
    % convert metadata into matlab struct from string
    % end of useful metadata
    endInd = strfind(metadata, '"RoiGroups"') - 5;
    % convert metadata to SI struct
    evalc(metadata(1:endInd));

    % get descriptions
    descriptions = imgData.descriptions();

    % save metadata and descriptions
    save('imMetaDat.mat', 'SI', 'descriptions', '-v7.3');

    % align time series

    % Create reference frame for image alignment 
    % reference stack using first numRefFrames frames 
    refStack = unalignedSeries(:,:,1:numRefFrames); 
    % max intensity projection of ref stack.
    refFrame = max(refStack, [], 3); 

    % Align the time series
    fprintf('Aligning %s \n', tifFileName);

    alignedSeries = fccAlignment(unalignedSeries, ...
        refFrame, 'xml');
    % Average over time of the whole aligned time series
    meanImageAligned = int16(mean(alignedSeries, 3));

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
    
    frameStartTimes = daqTime(frameStarts + 1);
    frameEndTimes = daqTime(frameEnds + 1);

    % save time series data
    save('imDat.mat', 'unalignedSeries', 'alignedSeries', ...
        'meanImageAligned', 'trialPath', 'tifFileName', ...
        'frameStartTimes', 'frameEndTimes', 'frameTimingError', '-v7.3');
                
end

% preprocessLegVid.m
%
% Function for preprocessing raw voltage signals from the leg tracking
%  camera (strobe signal with falling edge as frame start). Extracts leg
%  video frame times. Checks that number of frames camera says it has 
%  captured matches number of trigger pulses sent. Also, saves full path 
%  to leg tracking video.
%
%
% INPUTS:
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqOutput - struct of output signals sent by experimental DAQ,
%       processed by preprocessUserDaq() 
%   daqTime - vector of times corresponding to each sample point of daqData
%
% OUTPUTS:
%   None, but saves data in legVidDat.mat in pwd.
%
% CREATED: 12/4/18 HHY
% UPDATED: 2/25/19 HHY
%

function preprocessLegVid(daqData, daqOutput, daqTime)

    % frame start indicies, strobe signal - find falling edges
    frameStarts = find(diff(daqData.legCamFrames) < -0.1);
    
    % frame trigger indicies, output sent by experimental DAQ
    frameTrigs = find(diff(daqOutput.legCamFrameStartTrig) > 0.1);
    
    % check that number of captured frames matches number of triggered
    %  frames
    if (length(frameStarts) ~= length(frameTrigs))
        disp('Warning: Frame count mismatch in leg tracking video');
    end
    
    % leg vid frame times
    legVidFrameTimes = daqTime(frameStarts + 1);
    
    % get full path and file name of legVid - assumes that name contains
    %  legVid and the file is .mp4
    vidDir = dir('legVid*.mp4');
    
    % check that legVid video exists
    if (~isempty(vidDir))
        vidName = vidDir.name;
        vidFolder = vidDir.folder;
        vidExists = true;
    else
        disp('Warning: No leg video matching name legVid*.mp4 found');
        vidName = [];
        vidFolder = [];
        vidExists = false;
    end
    
    % save leg vid info in legVidDat.mat file in pwd
    save('legVidDat.mat', 'legVidFrameTimes', 'vidName', 'vidFolder',...
        'vidExists', '-v7.3');
    
    disp('Leg video data saved!');
    
end
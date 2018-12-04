% convertImgSeqToVideo.m
%
% Function that takes in path to a folder containing a sequence of images
%  and converts it to video of the specified type saved one level up from
%  the image sequence.
%
% NOTE: currently doesn't work correctly on default output of Point Gray
%  Spinnaker b/c file numbers aren't zero padded
%
% INPUTS:
%   imgSeqPath - full path to folder containing all the images
%   videoType - type of video, of profile options of VideoWriter
%   numVideos - number of videos to convert image sequence into
%   frameRate - frame rate of videos, in Hz
%
% OUTPUTS:
%   none, but generates videos saved one level up from image sequence as
%   side effect
%
% CREATED: 11/22/18
% UPDATED: 11/22/18 HHY
%

function convertImgSeqToVideo(imgSeqPath, videoType, numVideos, frameRate)

    curDir = pwd; % save present working directory
    
    % change to folder containing imgSeq folder
    cd(imgSeqPath)
    cd .. 

    % name of image sequence folder
    fsInd = strfind(imgSeqPath, filesep);
    imgSeqName = imgSeqPath((fsInd(end)+1):end);
    
    % get image files in folder
    imgSeq = dir(imgSeqPath);
    
    % 1st 2 of imgSeq are not images
    numImgInVid = ceil((length(imgSeq) - 2) / numVideos);
    
    % counter of which image we're on
    counter = 3; % start at 3 b/c 1st 2 of imgSeq are not images
    
    % loop through all videos to create
    for i = 1:numVideos
        % name of this video
        videoName = sprintf('%s_%02d', imgSeqName, i);
        
        % start VideoWriter
        v = VideoWriter(videoName, videoType);
        v.FrameRate = frameRate;
        open(v)
        
        disp(['Writing ' videoName]);
        % write individual images to video
        for j = 1:numImgInVid
            if (counter <= length(imgSeq))
                img = imread(...
                    [imgSeq(counter).folder filesep imgSeq(counter).name]);
                writeVideo(v,img)
                counter = counter + 1;
            end
        end
        
        close(v)
    end
end
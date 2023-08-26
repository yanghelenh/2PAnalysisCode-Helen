% getDFFfromBouts.m
%
% Helper function for saveDFFCond_bouts() that takes in dF/F and bout times
%  and returns the dF/F values aligned to the velocity peak
% To allow averaging over different bouts and the inconsistent frame rate,
%  interpolate to specified frame rate
%
% INPUTS:
%   dff - struct of imaging data, from img struct, with fields left, right,
%     sum, diff
%   imgT - time points for imaging data
%   dffParams - struct of parameters, directly from
%     saveDFFCond_bouts()
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to, in Hz
%   peakTime - start times of bout peaks
%   boutStartTime - start times of bout starts
%   boutEndTime - end times of bout ends
%
% OUTPUTS:
%   algnDFF - struct with fields left, right, sum, diff for all dF/F data,
%       each field is numTimePts x numBouts matrix for spikerate for bouts, 
%       aligned to yaw velocity peak
%   t - time vector for dF/F, interpolated
%
% CREATED: 8/25/23 - HHY
%
% UPDATED:
%   8/25/23 - HHY
%
function [algnDFF, t] = getDFFfromBouts(dff, imgT, ...
    dffParams, peakTime, boutStartTime, boutEndTime)

    % get number time pts to either side of peak, after interpolation
    maxNumFrames = floor(dffParams.maxDuration * ...
        dffParams.interpFrameRate);
    % number of bouts
    numBouts = length(peakTime);
    % interframe interval, for interpolation
    ifi = 1/dffParams.interpFrameRate;

    % preallocate
    algnDFF.left = nan(maxNumFrames * 2 + 1,numBouts);
    algnDFF.right = nan(maxNumFrames * 2 + 1,numBouts);
    algnDFF.sum = nan(maxNumFrames * 2 + 1,numBouts);
    algnDFF.diff = nan(maxNumFrames * 2 + 1,numBouts);


    % loop through all bouts
    for i = 1:numBouts
        % zero time vector, so that yaw velocity peaks are at t = 0
        % account for delay
        tOrig = imgT - peakTime(i);

        % get time for bout start and end, zeroed
        boutStartT = boutStartTime(i) - peakTime(i);
        boutEndT = boutEndTime(i) - peakTime(i);

        % get time vector for interpolation - only consider +/- maxDuration
        %  around peak
        % doing it this way keeps 0 at 0
        newTDur = (maxNumFrames / dffParams.interpFrameRate);


        newTHalf1 = 0:ifi:newTDur;
        newTHalf1 = fliplr(newTHalf1) * -1;

        newTHalf2 = 0:ifi:newTDur;

        newT = [newTHalf1 newTHalf2(2:end)];

        % get interpolated dF/F
        interpLeft = interp1(tOrig, dff.left, newT,'linear');
        interpRight = interp1(tOrig, dff.right, newT,'linear');
        interpDiff = interp1(tOrig, dff.diff, newT,'linear');
        interpSum = interp1(tOrig, dff.sum, newT,'linear');


        % add to output matrix
        algnDFF.left(:,i) = interpLeft';
        algnDFF.right(:,i) = interpRight';
        algnDFF.sum(:,i) = interpSum';
        algnDFF.diff(:,i) = interpDiff';
    end

    % get time vector for output matrix
    tHalf1 = fliplr((0:maxNumFrames) * ifi * -1);
    tHalf2 = (1:maxNumFrames) * ifi;

    t = [tHalf1 tHalf2]';    
end
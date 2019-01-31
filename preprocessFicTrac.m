% preprocessFicTrac.m
%
% Function for preprocessing raw FicTrac voltage data acquired on the user
%  DAQ. Called by preprocess(). Takes 0-10V signal of heading, intX, and
%  intY and returns angular position and velocity for each dimension. For x
%  and y, also returns position and velocity in distance units (mm). Raw
%  voltage signal is lowpass filtered 2X with cutoff specified by
%  LOWPASS_FILTER_CUTOFF. MAX_YAW_VELOCITY, MAX_FWD_VELOCITY, 
%  MAX_SLIDE_VELOCITY specify max anglular velocities in each of these 
%  dimensions. If they are exceeded, program assumes they're errors and 
%  replaces them with previous value that was below threshold.
%
% NOTE: Current output of FicTrac (using Wilson lab shenanigans version of
%  FicTrac.cpp), intX and intY are forward and slide (i.e. they don't
%  incorporate fly's heading)
%
%
% INPUTS:
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqTime - vector of times corresponding to each sample point of daqData
%   sampRate - sampling rate of acquisition
%
% OUTPUTS:
%   None, but saves data in fictracDat.mat in pwd.
%
% CREATED: 12/3/18 HHY
% UPDATED: 12/4/18 HHY
%

function preprocessFicTrac(daqData, daqTime, sampRate)
    % constants
    % lowpass filter cutoff (approximately half of FicTrac sample rate)
    LOWPASS_FILTER_CUTOFF = 40; % in Hz 
    MAX_YAW_VELOCITY = 2500; % deg/sec
    MAX_FWD_VELOCITY = 2500; % deg/sec
    MAX_SLIDE_VELOCITY = 2500; % deg/sec
    BALL_DIAM = 6.46; % diameter of ball, in mm
    
    % yaw/heading
    [yawAngVel, yawAngPos] = ficTracSignalDecoding(...
        daqData.ficTracHeading, sampRate, LOWPASS_FILTER_CUTOFF, ...
        MAX_YAW_VELOCITY);
    
    % wrap yaw angular position to 360 deg instead of it being cumulative
    yawAngPosWrap = wrapTo360(yawAngPos);
    
    % conversion factor between degrees and mm
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    mmPerDeg = circum / 360; % mm per degree of ball
    
    % forward direction (intX)
    [fwdAngVel, fwdAngPos] = ficTracSignalDecoding(daqData.ficTracIntX, ...
        sampRate, LOWPASS_FILTER_CUTOFF, MAX_FWD_VELOCITY);
    fwdVel = fwdAngVel .* mmPerDeg; % velocity in mm/sec
    % cumulative forward position in mm, where start of trial is at 0
    fwdCumPos = (fwdAngPos - fwdAngPos(1)) .* mmPerDeg; 
    
    % slide direction (intY)
    [slideAngVel, slideAngPos] = ficTracSignalDecoding(...
        daqData.ficTracIntY, sampRate, LOWPASS_FILTER_CUTOFF, ...
        MAX_SLIDE_VELOCITY);
    slideVel = slideAngVel .* mmPerDeg; % velocity in mm/sec
    % cumulative slide position in mm, where start of trial is at 0
    slideCumPos = (slideAngPos-slideAngPos(1)) .* mmPerDeg;     
    
    % position incorporating heading - as if fly were walking on x-y plane,
    %  x-y coordinates at each time point
    % start with fly at (0,0) and facing 0 deg
    zeroedYawAngPos = yawAngPos - yawAngPos(1); 
    
    % movement in x (in degrees) at each time point
    xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
        (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

    % x position in mm (i.e. x-coordinate of fly's position at each time 
    %  point), starts at 0
    xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
   
    % movement in y (in degrees) at each time point
    yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
        (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

    % y position in mm (i.e. y-coordinate of fly's position at each time 
    %  point), starts at 0
    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
    
    % time
    t = daqTime';
    
    % save data
    save('fictracDat.mat', 'yawAngVel', 'yawAngPosWrap', 'fwdVel', ...
        'fwdCumPos', 'slideVel', 'slideCumPos', 'xPos', 'yPos', ...
        't', '-v7.3');

end
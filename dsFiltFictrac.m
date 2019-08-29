% dsFiltFictrac.m
%
% Function that downsamples and smooths fictrac data and returns new fwd,
%  slide, and yaw postion and velocity, x and y position, total speed, and
%  yaw speed, as well as t and dropInd in fictrac struct
%
% INPUT:
%   trialPath - full path to trial folder that contains fictracDat.mat
%   fictrac - struct of fictrac parameters/data, when passed as input,
%       contains:
%           dsf - downsampling factor
%           filtParams - struct of parameters for smoothing all fictrac 
%                   data
%               padLen - padding length for convolving with Gaussian kernel
%               sigmaPos - standard deviation for Gaussian kernel for 
%                   position
%               sigmaVel - standard deviation for Gaussian kernel for 
%                   velocity
%           degPerMM - degrees per millimeter, conversion factor
%           mmPerDeg - millimeters per degree, conversion factor
%
% OUTPUT:
%   fictrac - struct of fictrac parameters/data, now also contains:
%       fwdCumPos
%       fwdVel
%       slideCumPos
%       slideVel
%       yawAngCumPos
%       yawAngPosWrap
%       yawAngVel
%       yawAngSpd
%       totAngSpd
%       totSpd
%       xPos
%       yPos
%       t
%       dropInd
%
% CREATED: 8/27/19 - HHY
%
% UPDATED: 8/27/19 - HHY
%

function fictrac = dsFiltFictrac(trialPath, fictrac)
    curDir = pwd;
    
    cd(trialPath);
    
    % load data
    load('fictracDat.mat', 'dropInd', 'fwdCumPos', 'slideCumPos', 't', ...
        'yawAngPosWrap');
    
    % downsample to 500 Hz, all position variables are cumulative
    fwdPosDS = downsample(fwdCumPos, fictrac.dsf);
    yawAngPos = unwrap(yawAngPosWrap .* (pi / 180));
    yawAngPosDS = downsample(yawAngPos, fictrac.dsf);
    yawAngPosDS = yawAngPosDS .* (180 / pi);
    slidePosDS = downsample(slideCumPos, fictrac.dsf);
    
    % downsample timing vector and dropInd to 500 Hz
    tDS = downsample(t, fictrac.dsf);
    ftSampRate = 1/median(diff(tDS));

    dropIndDS = unique(round(dropInd ./ fictrac.dsf));
    % remove any zeros from dropIndDS
    dropIndDS(dropIndDS < 1) = [];
    
    % perform smoothing on position and extract smoothed velocity
    [fwdPosSmo, fwdVelSmo] = computeSmoothedVelocity(fwdPosDS, ...
        fictrac.filtParams.padLen, fictrac.filtParams.sigmaPos, ...
        fictrac.filtParams.sigmaVel, ftSampRate);
    [yawAngPosSmo, yawAngVelSmo] = computeSmoothedVelocity(yawAngPosDS, ...
        fictrac.filtParams.padLen, fictrac.filtParams.sigmaPos, ...
        fictrac.filtParams.sigmaVel, ftSampRate);
    [slidePosSmo, slideVelSmo] = computeSmoothedVelocity(slidePosDS, ...
        fictrac.filtParams.padLen, fictrac.filtParams.sigmaPos, ...
        fictrac.filtParams.sigmaVel, ftSampRate);
    
    % compute wrapped yaw position
    yawAngPosSmoWrap = wrapTo360(yawAngPosSmo);
    
    % compute yaw speed
    yawAngSpd = abs(yawAngVelSmo);
    
    % compute total speed, in degrees per second
    % convert forward and slide velocities from mm to deg per sec
    fwdVelDeg = fwdVelSmo .* fictrac.degPerMM;
    slideVelDeg = slideVelSmo .* fictrac.degPerMM;
    
    totAngSpd = yawAngSpd + abs(fwdVelDeg) + abs(slideVelDeg);
    totSpd = totAngSpd .* fictrac.mmPerDeg;
    
    % compute x and y positions, using smoothed positions and velocities
    % start with fly at (0,0) and facing 0 deg
    zeroedYawAngPos = yawAngPosSmo - yawAngPosSmo(1); 

    % movement in x (in degrees) at each time point
    xChangePos = (fwdVelDeg ./ ftSampRate) .* sind(zeroedYawAngPos) + ...
        (slideVelDeg ./ ftSampRate) .* sind(zeroedYawAngPos + 90);  

    % x position in mm (i.e. x-coordinate of fly's position at each time 
    %  point), starts at 0
    xPos = (cumsum(xChangePos) - xChangePos(1)) .* fictrac.mmPerDeg;
   
    % movement in y (in degrees) at each time point
    yChangePos = (fwdVelDeg ./ ftSampRate) .* cosd(zeroedYawAngPos) + ...
        (slideVelDeg ./ ftSampRate) .* cosd(zeroedYawAngPos + 90);

    % y position in mm (i.e. y-coordinate of fly's position at each time 
    %  point), starts at 0
    yPos = (cumsum(yChangePos) - yChangePos(1)) .* fictrac.mmPerDeg;
    
    
    % save all data into fictrac struct
    fictrac.fwdCumPos = fwdPosSmo;
    fictrac.fwdVel = fwdVelSmo;
    fictrac.slideCumPos = slidePosSmo;
    fictrac.slideVel = slideVelSmo;
    fictrac.yawAngCumPos = yawAngPosSmo;
    fictrac.yawAngPosWrap = yawAngPosSmoWrap;
    fictrac.yawAngVel = yawAngVelSmo;
    fictrac.yawAngSpd = yawAngSpd;
    fictrac.totAngSpd = totAngSpd;
    fictrac.totSpd = totSpd;
    fictrac.xPos = xPos;
    fictrac.yPos = yPos;
    fictrac.t = tDS'; % change to row vector
    fictrac.dropInd = dropIndDS;
    
    cd(curDir);

end
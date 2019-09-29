% addAcceleration2PData.m
%
% Function that loops through all pData in folder provided, computes
%  acceleration from FicTrac velocity measurements, and saves it back into
%  the pData.
% Smooths velocity again before computing acceleration and smooths
%  acceleration before saving it.
%
% INPUTS:
%   pDat - full path to pData directory
%   sigmaVel2 - standard deviation for Gaussian kernel for smoothing for
%       velocity again
%   sigmaAcc - standard deviation for Gaussian kernel for smoothing for
%       acceleration
%
% OUTPUTS: none, but resaves pData file with updated fictrac struct:
%   filtParams.sigmaVel2 - from input
%   filtParams.sigmaAcc - from input
%   fwdAcc - forward acceleration (in mm/s2)
%   yawAngAcc - yaw acceleration (in deg/s2)
%   slideAcc - slide acceleration (in mm/s2)
%   totAngAccMag - total magnitude of acceleration in all three axes (in
%       deg/s2)
%   totAccMag - same as above but in mm/s2
%
% CREATED: 9/27/19 - HHY
%
% UPDATED:
%   9/27/19 - HHY
%

function addAcceleration2PData(pDat, sigmaVel2, sigmaAcc)

    pDataFiles = dir([pDat filesep '*_pData.mat']);
    
    for i = 1:length(pDataFiles)
        pDatFullPath = [pDataFiles(i).folder filesep pDataFiles(i).name];
        
        load(pDatFullPath, 'fictrac');
        
        sampRate = 1/median(diff(fictrac.t));
        
        % compute acceleration (function works to smooth and take the
        %  derivative, so will return acceleration)
        [~, fictrac.fwdAcc] = computeSmoothedVelocity(fictrac.fwdVel, ...
            fictrac.filtParams.padLen, sigmaVel2, sigmaAcc, sampRate);
        [~, fictrac.yawAngAcc] = computeSmoothedVelocity(fictrac.yawAngVel, ...
            fictrac.filtParams.padLen, sigmaVel2, sigmaAcc, sampRate);
        [~, fictrac.slideAcc] = computeSmoothedVelocity(fictrac.slideVel, ...
            fictrac.filtParams.padLen, sigmaVel2, sigmaAcc, sampRate);
        
        % convert to deg
        fwdAccDeg = fictrac.fwdAcc .* fictrac.degPerMM;
        slideAccDeg = fictrac.slideAcc .* fictrac.degPerMM;
        
        % compute magnitude of acceleration
        fictrac.totAngAccMag = abs(fwdAccDeg) + abs(slideAccDeg) + ...
            abs(fictrac.yawAngAcc);
        fictrac.totAccMag = fictrac.totAngAccMag .* fictrac.mmPerDeg;
        
        % save back into pData file
        save(pDatFullPath, 'fictrac', '-append');
    end
end
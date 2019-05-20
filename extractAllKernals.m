% extractAllKernals.m
%
% Function to calculate, plot, and save kernals relating responses of ROIs
%  to FicTrac variables.
%
% Make sure to have calculated dF/F and identified FicTrac dropped frames
%  before running this function
%
% INPUT:
%   winLen - length of window, in seconds that averages are computed over
%   lowPassCutoff - cutoff of lowpass filter to apply to numerator, use 0 
%       if no filtering desired
%   sampRate - sampling rate to convert dF/F and FicTrac data to, and to
%       calculate kernal at
%   also prompts user for trial folder
%
% OUTPUT:
%   none - but saves kernals to pData file and generates plots of kernals
%
% CREATED: 4/15/19 HHY
%
% UPDATED: 4/15/19 HHY
%

function extractAllKernals(winLen, lowpassCutoff, sampRate)

    % ask user to select trial folder
    disp('Select a trial folder to display.');
    uTrialPath = uigetdir;
    curDir = pwd;
    cd(uTrialPath)
    
    fprintf('Displaying %s \n', uTrialPath);
    
    % load relevant data
    load('pData.mat', 'dFFs', 'frameStartTimes');
    load('fictracDat.mat', 'dropInd', 'fwdVel', 't', 'yawAngVel');
    
    % turn FicTrac values to NaNs where it dropped
    fwdVel(dropInd) = nan; 
    yawAngVel(dropInd) = nan;
    
    % get yaw speed
    yawAngSpd = abs(yawAngVel);
    
    % number of ROIs
    numROIs = size(dFFs,1);
    
    % get sum and difference of dF/F responses; only for 2 ROIs
    if (numROIs == 2)
        sumDFF = dFFs(1,:) + dFFs(2,:);
        diffDFF = dFFs(2,:) - dFFs(1,:); % right - left ROIs
    end
    
    % convert dF/F data and FicTrac data to same timescale
    newT = frameStartTimes(1):(1/sampRate):frameStartTimes(end);
    % linear interpolation of data
    dFFsRS = interp1(frameStartTimes, dFFs', newT)';
    fwdVelRS = interp1(t, fwdVel, newT);
    yawAngVelRS = interp1(t, yawAngVel, newT);
    yawAngSpdRS = interp1(t, yawAngSpd, newT);
    sumDFFrs = interp1(frameStartTimes, sumDFF, newT);
    diffDFFrs = interp1(frameStartTimes, diffDFF, newT);
    
    % convert window from seconds to samples
    windowSamps = winLen * sampRate;
    
       
    % compute kernals for all ROIs
    for i = 1:numROIs
        % forward kernal - fwd vel and dF/F
        kernalsIndiv(i).fFwdVel = computeWienerKernel(fwdVelRS, dFFsRS(i,:),...
            sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernal - fwd vel and dF/F
        kernalsIndiv(i).rFwdVel = computeWienerKernel(dFFsRS(i,:), fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff);
        
        % forward kernal - yaw vel and dF/F
        kernalsIndiv(i).fYawVel = computeWienerKernel(yawAngVelRS, ...
            dFFsRS(i,:), sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernal - yaw vel and dF/F
        kernalsIndiv(i).rYawVel = computeWienerKernel(dFFsRS(i,:), ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        
        % forward kernal - yaw speed and dF/F
        kernalsIndiv(i).fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            dFFsRS(i,:), sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernal - yaw speed and dF/F
        kernalsIndiv(i).rYawSpd = computeWienerKernel(dFFsRS(i,:), ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);
    end
    
    % compute sum and diff kernals
    if (numROIs == 2)
        % sum kernals
        kernalsSum.fFwdVel = computeWienerKernel(fwdVelRS, sumDFFrs,...
            sampRate, windowSamps, lowpassCutoff);
        kernalsSum.rFwdVel = computeWienerKernel(sumDFFrs, fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff); 
        kernalsSum.fYawVel = computeWienerKernel(yawAngVelRS, ...
            sumDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernalsSum.rYawVel = computeWienerKernel(sumDFFrs, ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        kernalsSum.fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            sumDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernalsSum.rYawSpd = computeWienerKernel(sumDFFrs, ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);
        
        % diff kernals
        kernalsDiff.fFwdVel = computeWienerKernel(fwdVelRS, diffDFFrs,...
            sampRate, windowSamps, lowpassCutoff);
        kernalsDiff.rFwdVel = computeWienerKernel(diffDFFrs, fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff); 
        kernalsDiff.fYawVel = computeWienerKernel(yawAngVelRS, ...
            diffDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernalsDiff.rYawVel = computeWienerKernel(diffDFFrs, ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        kernalsDiff.fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            diffDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernalsDiff.rYawSpd = computeWienerKernel(diffDFFrs, ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);        
        
    end
    
    % kernal timescale
    kT = 0:1:(windowSamps - 1);
    kT = kT * (1/sampRate);
    
    % generate plots
    for i = 1:numROIs
        figure;
        subplot(2, 3, 1);
        plot(kT, kernalsIndiv(i).fFwdVel);
        title('Fwd Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernalsIndiv(i).fYawVel);
        title('Fwd Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernalsIndiv(i).fYawSpd);
        title('Fwd Kernal - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernalsIndiv(i).rFwdVel);
        title('Rev Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernalsIndiv(i).rYawVel);
        title('Rev Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernalsIndiv(i).rYawSpd);
        title('Rev Kernal - Yaw Speed');
        xlabel('time (s)');
        
    end
    
    if (numROIs == 2)
        figure; 
        
        subplot(2, 3, 1);
        plot(kT, kernalsSum.fFwdVel);
        title('Fwd Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernalsSum.fYawVel);
        title('Fwd Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernalsSum.fYawSpd);
        title('Fwd Kernal - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernalsSum.rFwdVel);
        title('Rev Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernalsSum.rYawVel);
        title('Rev Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernalsSum.rYawSpd);
        title('Rev Kernal - Yaw Speed');
        xlabel('time (s)');
        
        figure; 
        
        subplot(2, 3, 1);
        plot(kT, kernalsDiff.fFwdVel);
        title('Fwd Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernalsDiff.fYawVel);
        title('Fwd Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernalsDiff.fYawSpd);
        title('Fwd Kernal - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernalsDiff.rFwdVel);
        title('Rev Kernal - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernalsDiff.rYawVel);
        title('Rev Kernal - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernalsDiff.rYawSpd);
        title('Rev Kernal - Yaw Speed');
        xlabel('time (s)');        
    end
    
    % save kernal parameters in struct
    kernalParams.t = kT;
    kernalParams.winLen = winLen;
    kernalParams.lowpassCutoff = lowpassCutoff;
    kernalParams.sampRate = sampRate;
    
    % save kernals into pData file
    if (numROIs == 2)
        save('pData.mat', 'kernalsIndiv', 'kernalsSum', 'kernalsDiff', ...
        	'kernalParams', '-append');
    else
        save('pData.mat', 'kernalsIndiv', 'kernalParams', '-append');    
    end  
end
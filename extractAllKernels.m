% extractAllkernels.m
%
% Function to calculate, plot, and save kernels relating responses of ROIs
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
%       calculate kernel at
%   also prompts user for trial folder
%
% OUTPUT:
%   none - but saves kernels to pData file and generates plots of kernels
%
% CREATED: 4/15/19 HHY
%
% UPDATED: 4/15/19 HHY
%

function extractAllKernels(winLen, lowpassCutoff, sampRate)

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
    
    if (numROIs == 2)
        sumDFFrs = interp1(frameStartTimes, sumDFF, newT);
        diffDFFrs = interp1(frameStartTimes, diffDFF, newT);
    elseif (numROIs == 1)
        dFFsRS = dFFsRS'; % otherwise, indexing wrong
    end
    
    % convert window from seconds to samples
    windowSamps = winLen * sampRate;
    
       
    % compute kernels for all ROIs
    for i = 1:numROIs
        % forward kernel - fwd vel and dF/F
        kernelsIndiv(i).fFwdVel = computeWienerKernel(fwdVelRS, dFFsRS(i,:),...
            sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernel - fwd vel and dF/F
        kernelsIndiv(i).rFwdVel = computeWienerKernel(dFFsRS(i,:), fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff);
        
        % forward kernel - yaw vel and dF/F
        kernelsIndiv(i).fYawVel = computeWienerKernel(yawAngVelRS, ...
            dFFsRS(i,:), sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernel - yaw vel and dF/F
        kernelsIndiv(i).rYawVel = computeWienerKernel(dFFsRS(i,:), ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        
        % forward kernel - yaw speed and dF/F
        kernelsIndiv(i).fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            dFFsRS(i,:), sampRate, windowSamps, lowpassCutoff);
        
        % reverse kernel - yaw speed and dF/F
        kernelsIndiv(i).rYawSpd = computeWienerKernel(dFFsRS(i,:), ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);
    end
    
    % compute sum and diff kernels
    if (numROIs == 2)
        % sum kernels
        kernelsSum.fFwdVel = computeWienerKernel(fwdVelRS, sumDFFrs,...
            sampRate, windowSamps, lowpassCutoff);
        kernelsSum.rFwdVel = computeWienerKernel(sumDFFrs, fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff); 
        kernelsSum.fYawVel = computeWienerKernel(yawAngVelRS, ...
            sumDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernelsSum.rYawVel = computeWienerKernel(sumDFFrs, ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        kernelsSum.fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            sumDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernelsSum.rYawSpd = computeWienerKernel(sumDFFrs, ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);
        
        % diff kernels
        kernelsDiff.fFwdVel = computeWienerKernel(fwdVelRS, diffDFFrs,...
            sampRate, windowSamps, lowpassCutoff);
        kernelsDiff.rFwdVel = computeWienerKernel(diffDFFrs, fwdVelRS,...
            sampRate, windowSamps, lowpassCutoff); 
        kernelsDiff.fYawVel = computeWienerKernel(yawAngVelRS, ...
            diffDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernelsDiff.rYawVel = computeWienerKernel(diffDFFrs, ...
            yawAngVelRS, sampRate, windowSamps, lowpassCutoff);
        kernelsDiff.fYawSpd = computeWienerKernel(yawAngSpdRS, ...
            diffDFFrs, sampRate, windowSamps, lowpassCutoff);
        kernelsDiff.rYawSpd = computeWienerKernel(diffDFFrs, ...
            yawAngSpdRS, sampRate, windowSamps, lowpassCutoff);        
        
    end
    
    % kernel timescale
    kT = 0:1:(windowSamps - 1);
    kT = kT * (1/sampRate);
    
    % generate plots
    for i = 1:numROIs
        figure;
        subplot(2, 3, 1);
        plot(kT, kernelsIndiv(i).fFwdVel);
        title('Fwd kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernelsIndiv(i).fYawVel);
        title('Fwd kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernelsIndiv(i).fYawSpd);
        title('Fwd kernel - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernelsIndiv(i).rFwdVel);
        title('Rev kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernelsIndiv(i).rYawVel);
        title('Rev kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernelsIndiv(i).rYawSpd);
        title('Rev kernel - Yaw Speed');
        xlabel('time (s)');
        
    end
    
    if (numROIs == 2)
        figure; 
        
        subplot(2, 3, 1);
        plot(kT, kernelsSum.fFwdVel);
        title('Fwd kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernelsSum.fYawVel);
        title('Fwd kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernelsSum.fYawSpd);
        title('Fwd kernel - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernelsSum.rFwdVel);
        title('Rev kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernelsSum.rYawVel);
        title('Rev kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernelsSum.rYawSpd);
        title('Rev kernel - Yaw Speed');
        xlabel('time (s)');
        
        figure; 
        
        subplot(2, 3, 1);
        plot(kT, kernelsDiff.fFwdVel);
        title('Fwd kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 2);
        plot(kT, kernelsDiff.fYawVel);
        title('Fwd kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 3);
        plot(kT, kernelsDiff.fYawSpd);
        title('Fwd kernel - Yaw Speed');
        xlabel('time (s)');
        
        subplot(2, 3, 4);
        plot(kT, kernelsDiff.rFwdVel);
        title('Rev kernel - Fwd Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 5);
        plot(kT, kernelsDiff.rYawVel);
        title('Rev kernel - Yaw Velocity');
        xlabel('time (s)');
        
        subplot(2, 3, 6);
        plot(kT, kernelsDiff.rYawSpd);
        title('Rev kernel - Yaw Speed');
        xlabel('time (s)');        
    end
    
    % save kernel parameters in struct
    kernelParams.t = kT;
    kernelParams.winLen = winLen;
    kernelParams.lowpassCutoff = lowpassCutoff;
    kernelParams.sampRate = sampRate;
    
    % save kernels into pData file
    if (numROIs == 2)
        save('pData.mat', 'kernelsIndiv', 'kernelsSum', 'kernelsDiff', ...
        	'kernelParams', '-append');
        disp('Saved!');
    else
        save('pData.mat', 'kernelsIndiv', 'kernelParams', '-append');
        disp('Saved!');
    end  
    cd(curDir);
end
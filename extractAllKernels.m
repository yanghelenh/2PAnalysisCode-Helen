% extractAllkernels.m
%
% Function to calculate, plot, and save kernels relating responses of ROIs
%  to FicTrac variables. Calls computeWienerKernel()
%
% Make sure to have calculated dF/F and identified FicTrac dropped frames
%  before running this function
%
% INPUT:
%   trialPath - full path to trial folder to extract kernel
%   winLen - length of window, in seconds that averages are computed over
%   cutFreq - cutoff frequency (f_cut) in attenuation applied to frequency
%       domain filter, a la Nagel and Wilson 2011. Set to 0 with tauFreq if
%       no attenuation desired.
%   tauFreq - f_tau in attenuation applied to frequency domain filter, a la
%       Nagel and Wilson 2011. Set to 0 with cutFreq if no attenuation
%       desired.
%   sampRate - sampling rate to convert dF/F and FicTrac data to, and to
%       calculate kernel at
%   doDisp - binary flag for whether or not to plot kernels
%
% OUTPUT:
%   none - but saves kernels to pData file and generates plots of kernels
%       if doDisp == 1
%
% CREATED: 4/15/19 HHY
% UPDATED: 5/21/19 HHY
%

function extractAllKernels(trialPath, winLen, cutFreq, tauFreq, ...
    sampRate, doDisp)
    
    curDir = pwd;
    
    cd(trialPath)
    
    trialPathFiles = dir(trialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasPDat = sum(strcmp(trialPathFileNames, 'pData.mat'));
    hasFTDat = sum(strcmp(trialPathFileNames, 'fictracDat.mat'));
    
    if (hasPDat && hasFTDat)
        % load relevant data
        load('pData.mat', 'dFFs', 'frameStartTimes');
        load('fictracDat.mat', 'dropInd', 'fwdVel', 't', 'yawAngVel');
        
        % temporarily test timing with frameEndTimes
%         load('imDat.mat','frameEndTimes');
%         frameStartTimes = frameEndTimes;

        % turn FicTrac values to NaNs where it dropped
        fwdVel(dropInd) = nan; 
        yawAngVel(dropInd) = nan;

        % get yaw speed
        yawAngSpd = abs(yawAngVel);

        % number of ROIs
        numROIs = size(dFFs,1);

        % get sum and difference of dF/F responses; only for 2 ROIs
        % z-score each so left and right cells are normalized
        if (numROIs == 2)
%             zLeft = zscore(dFFs(1,:));
%             zRight = zscore(dFFs(2,:));
%             sumDFF = zLeft + zRight;
%             diffDFF = zRight - zLeft; % right - left ROIs
            sumDFF = dFFs(1,:) + dFFs(2,:);
            diffDFF = dFFs(2,:) - dFFs(1,:);
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
            [kernelsIndiv(i).fFwdVel, ~] = computeWienerKernel(...
                fwdVelRS, dFFsRS(i,:), sampRate, winLen, cutFreq, tauFreq);

            % reverse kernel - fwd vel and dF/F
            [kernelsIndiv(i).rFwdVel, ~] = computeWienerKernel(...
                dFFsRS(i,:),fwdVelRS, sampRate, winLen, cutFreq, tauFreq);

            % forward kernel - yaw vel and dF/F
            [kernelsIndiv(i).fYawVel,~] = computeWienerKernel(...
                yawAngVelRS, dFFsRS(i,:), sampRate, winLen, cutFreq,...
                tauFreq);

            % reverse kernel - yaw vel and dF/F
            [kernelsIndiv(i).rYawVel,~] = computeWienerKernel(...
                dFFsRS(i,:), yawAngVelRS, sampRate, winLen, cutFreq, ...
                tauFreq);

            % forward kernel - yaw speed and dF/F
            [kernelsIndiv(i).fYawSpd, ~] = computeWienerKernel(...
                yawAngSpdRS, dFFsRS(i,:), sampRate, winLen, cutFreq, ...
                tauFreq);

            % reverse kernel - yaw speed and dF/F
            [kernelsIndiv(i).rYawSpd, lags] = computeWienerKernel(...
                dFFsRS(i,:), yawAngSpdRS, sampRate, winLen, cutFreq, ...
                tauFreq);
        end

        % compute sum and diff kernels
        if (numROIs == 2)
            % sum kernels
            [kernelsSum.fFwdVel, ~] = computeWienerKernel(fwdVelRS, ...
                sumDFFrs,sampRate, winLen, cutFreq, tauFreq);
            [kernelsSum.rFwdVel, ~] = computeWienerKernel(sumDFFrs, ...
                fwdVelRS,sampRate, winLen, cutFreq, tauFreq); 
            [kernelsSum.fYawVel, ~] = computeWienerKernel(yawAngVelRS, ...
                sumDFFrs, sampRate, winLen, cutFreq, tauFreq);
            [kernelsSum.rYawVel, ~] = computeWienerKernel(sumDFFrs, ...
                yawAngVelRS, sampRate, winLen, cutFreq, tauFreq);
            [kernelsSum.fYawSpd, ~] = computeWienerKernel(yawAngSpdRS, ...
                sumDFFrs, sampRate, winLen, cutFreq, tauFreq);
            [kernelsSum.rYawSpd, ~] = computeWienerKernel(sumDFFrs, ...
                yawAngSpdRS, sampRate, winLen, cutFreq, tauFreq);

            % diff kernels
            [kernelsDiff.fFwdVel, ~] = computeWienerKernel(fwdVelRS, ...
                diffDFFrs, sampRate, winLen, cutFreq, tauFreq);
            [kernelsDiff.rFwdVel, ~] = computeWienerKernel(diffDFFrs, ...
                fwdVelRS, sampRate, winLen, cutFreq, tauFreq); 
            [kernelsDiff.fYawVel, ~] = computeWienerKernel(yawAngVelRS, ...
                diffDFFrs, sampRate, winLen, cutFreq, tauFreq);
            [kernelsDiff.rYawVel, ~] = computeWienerKernel(diffDFFrs, ...
                yawAngVelRS, sampRate, winLen, cutFreq, tauFreq);
            [kernelsDiff.fYawSpd, ~] = computeWienerKernel(yawAngSpdRS, ...
                diffDFFrs, sampRate, winLen, cutFreq, tauFreq);
            [kernelsDiff.rYawSpd, ~] = computeWienerKernel(diffDFFrs, ...
                yawAngSpdRS, sampRate, winLen, cutFreq, tauFreq);        

        end


        if doDisp
            % generate plots
            for i = 1:numROIs
                figure;
                subplot(2, 3, 1);
                plot(lags, kernelsIndiv(i).fFwdVel);
                title('Fwd kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 2);
                plot(lags, kernelsIndiv(i).fYawVel);
                title('Fwd kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 3);
                plot(lags, kernelsIndiv(i).fYawSpd);
                title('Fwd kernel - Yaw Speed');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 4);
                plot(lags, fliplr(kernelsIndiv(i).rFwdVel));
                title('Rev kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 5);
                plot(lags, fliplr(kernelsIndiv(i).rYawVel));
                title('Rev kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 6);
                plot(lags, fliplr(kernelsIndiv(i).rYawSpd));
                title('Rev kernel - Yaw Speed');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

            end

            if (numROIs == 2)
                figure; 

                subplot(2, 3, 1);
                plot(lags, kernelsSum.fFwdVel);
                title('Fwd kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 2);
                plot(lags, kernelsSum.fYawVel);
                title('Fwd kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 3);
                plot(lags, kernelsSum.fYawSpd);
                title('Fwd kernel - Yaw Speed');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 4);
                plot(lags, fliplr(kernelsSum.rFwdVel));
                title('Rev kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 5);
                plot(lags, fliplr(kernelsSum.rYawVel));
                title('Rev kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 6);
                plot(lags, fliplr(kernelsSum.rYawSpd));
                title('Rev kernel - Yaw Speed');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                figure; 

                subplot(2, 3, 1);
                plot(lags, kernelsDiff.fFwdVel);
                title('Fwd kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 2);
                plot(lags, kernelsDiff.fYawVel);
                title('Fwd kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 3);
                plot(lags, kernelsDiff.fYawSpd);
                title('Fwd kernel - Yaw Speed');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 4);
                plot(lags, fliplr(kernelsDiff.rFwdVel));
                title('Rev kernel - Fwd Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 5);
                plot(lags, fliplr(kernelsDiff.rYawVel));
                title('Rev kernel - Yaw Velocity');
                xlabel('time (s)');
                xlim([-1*winLen, winLen]);

                subplot(2, 3, 6);
                plot(lags, fliplr(kernelsDiff.rYawSpd));
                title('Rev kernel - Yaw Speed');
                xlabel('time (s)'); 
                xlim([-1*winLen, winLen]);
            end
        end

        % save kernel parameters in struct
        kernelParams.t = lags;
        kernelParams.winLen = winLen;
        kernelParams.cutFreq = cutFreq;
        kernelParams.tauFreq = tauFreq;
        kernelParams.sampRate = sampRate;

        % save kernels into pData file
    %     if (numROIs == 2)
    %         save('pData.mat', 'kernelsIndiv', 'kernelsSum', 'kernelsDiff', ...
    %         	'kernelParams', '-append');
    %         disp('Saved!');
    %     else
    %         save('pData.mat', 'kernelsIndiv', 'kernelParams', '-append');
    %         disp('Saved!');
    %     end  
        cd(curDir);
    else
        disp(['Selected trial folder does not contain pDat.mat file' ...
            ' and/or fictracDat.mat']);
        cd(curDir);
        return;
    end
end
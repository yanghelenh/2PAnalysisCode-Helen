% extractAllkernels.m
%
% Function to calculate, plot, and save kernels relating responses of ROIs
%  to FicTrac variables. Also, uses kernels to predict responses (in both
%  directions, Ca->FicTrac and FicTrac->Ca) and saves predicted responses
%  and variance explained.
% Calls computeWienerKernel(), computeLinearPrediction(),
%  computeVarianceExplained()
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
% UPDATED: 5/23/19 HHY
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
        fprintf('Kernels for \n%s \n', trialPath);
        
        % load relevant data
%         load('pData.mat', 'dFFs', 'frameStartTimes');
        
        % temporarily, to overwrite old kernels saved in pData; when done,
        %  replace with previous line
        load('pData.mat', 'dFFs', 'frameStartTimes','numChannels');
        
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
        if (numROIs == 2)
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


        % compute kernels, linear predictions, variance explained for 
        %  all ROIs
        for i = 1:numROIs
            % forward kernel - fwd vel and dF/F
            [kernels(i).fFwdVel.kernel, ~] = computeWienerKernel(...
                fwdVelRS, dFFsRS(i,:), sampRate, winLen, cutFreq, tauFreq);
            kernels(i).fFwdVel.predResp = computeLinearPrediction(...
                kernels(i).fFwdVel.kernel, fwdVelRS);
            [kernels(i).fFwdVel.varExpl, kernels(i).fFwdVel.fitObj] = ...
                computeVarianceExplained(dFFsRS(i,:), ...
                kernels(i).fFwdVel.predResp, 'poly1');

            % reverse kernel - fwd vel and dF/F
            [kernels(i).rFwdVel.kernel, ~] = computeWienerKernel(...
                dFFsRS(i,:), fwdVelRS, sampRate, winLen, cutFreq, tauFreq);
            kernels(i).rFwdVel.predResp = computeLinearPrediction(...
                kernels(i).rFwdVel.kernel, dFFsRS(i,:));
            [kernels(i).rFwdVel.varExpl, kernels(i).rFwdVel.fitObj] = ...
                computeVarianceExplained(fwdVelRS, ...
                kernels(i).rFwdVel.predResp, 'poly1');            

            % forward kernel - yaw vel and dF/F
            [kernels(i).fYawVel.kernel,~] = computeWienerKernel(...
                yawAngVelRS, dFFsRS(i,:), sampRate, winLen, cutFreq,...
                tauFreq);
            kernels(i).fYawVel.predResp = computeLinearPrediction(...
                kernels(i).fYawVel.kernel, yawAngVelRS);
            [kernels(i).fYawVel.varExpl, kernels(i).fYawVel.fitObj] = ...
                computeVarianceExplained(dFFsRS(i,:), ...
                kernels(i).fYawVel.predResp, 'poly1');            
            
            % reverse kernel - yaw vel and dF/F
            [kernels(i).rYawVel.kernel,~] = computeWienerKernel(...
                dFFsRS(i,:), yawAngVelRS, sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(i).rYawVel.predResp = computeLinearPrediction(...
                kernels(i).rYawVel.kernel, dFFsRS(i,:));
            [kernels(i).rYawVel.varExpl, kernels(i).rYawVel.fitObj] = ...
                computeVarianceExplained(yawAngVelRS, ...
                kernels(i).rYawVel.predResp, 'poly1'); 

            % forward kernel - yaw speed and dF/F
            [kernels(i).fYawSpd.kernel, ~] = computeWienerKernel(...
                yawAngSpdRS, dFFsRS(i,:), sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(i).fYawSpd.predResp = computeLinearPrediction(...
                kernels(i).fYawSpd.kernel, yawAngSpdRS);
            [kernels(i).fYawSpd.varExpl, kernels(i).fYawSpd.fitObj] = ...
                computeVarianceExplained(dFFsRS(i,:), ...
                kernels(i).fYawSpd.predResp, 'poly1');  

            % reverse kernel - yaw speed and dF/F
            [kernels(i).rYawSpd.kernel, lags] = computeWienerKernel(...
                dFFsRS(i,:), yawAngSpdRS, sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(i).rYawSpd.predResp = computeLinearPrediction(...
                kernels(i).rYawSpd.kernel, dFFsRS(i,:));
            [kernels(i).rYawSpd.varExpl, kernels(i).rYawSpd.fitObj] = ...
                computeVarianceExplained(yawAngSpdRS, ...
                kernels(i).rYawSpd.predResp, 'poly1'); 
            
        end

        % compute sum and diff kernels, linear predictions, 
        %  variance explained
        if (numROIs == 2)
            % sum kernels, predictions, variance; 3rd index
            [kernels(3).fFwdVel.kernel, ~] = computeWienerKernel(...
                fwdVelRS, sumDFFrs, sampRate, winLen, cutFreq, tauFreq);
            kernels(3).fFwdVel.predResp = computeLinearPrediction(...
                kernels(3).fFwdVel.kernel, fwdVelRS);
            [kernels(3).fFwdVel.varExpl, kernels(3).fFwdVel.fitObj] = ...
                computeVarianceExplained(sumDFFrs, ...
                kernels(3).fFwdVel.predResp, 'poly1');           
            
            [kernels(3).rFwdVel.kernel, ~] = computeWienerKernel(...
                sumDFFrs, fwdVelRS, sampRate, winLen, cutFreq, tauFreq); 
            kernels(3).rFwdVel.predResp = computeLinearPrediction(...
                kernels(3).rFwdVel.kernel, sumDFFrs);
            [kernels(3).rFwdVel.varExpl, kernels(3).rFwdVel.fitObj] = ...
                computeVarianceExplained(fwdVelRS, ...
                kernels(3).rFwdVel.predResp, 'poly1');               
            
            [kernels(3).fYawVel.kernel, ~] = computeWienerKernel(...
                yawAngVelRS, sumDFFrs, sampRate, winLen, cutFreq, tauFreq);
            kernels(3).fYawVel.predResp = computeLinearPrediction(...
                kernels(3).fYawVel.kernel, yawAngVelRS);
            [kernels(3).fYawVel.varExpl, kernels(3).fYawVel.fitObj] = ...
                computeVarianceExplained(sumDFFrs, ...
                kernels(3).fYawVel.predResp, 'poly1');              
            
            [kernels(3).rYawVel.kernel, ~] = computeWienerKernel(...
                sumDFFrs, yawAngVelRS, sampRate, winLen, cutFreq, tauFreq);
            kernels(3).rYawVel.predResp = computeLinearPrediction(...
                kernels(3).rYawVel.kernel, sumDFFrs);
            [kernels(3).rYawVel.varExpl, kernels(3).rYawVel.fitObj] = ...
                computeVarianceExplained(yawAngVelRS, ...
                kernels(3).rYawVel.predResp, 'poly1');                  
            
            [kernels(3).fYawSpd.kernel, ~] = computeWienerKernel(...
                yawAngSpdRS, sumDFFrs, sampRate, winLen, cutFreq, tauFreq);
            kernels(3).fYawSpd.predResp = computeLinearPrediction(...
                kernels(3).fYawSpd.kernel, yawAngSpdRS);
            [kernels(3).fYawSpd.varExpl, kernels(3).fYawSpd.fitObj] = ...
                computeVarianceExplained(sumDFFrs, ...
                kernels(3).fYawSpd.predResp, 'poly1');              
            
            [kernels(3).rYawSpd.kernel, ~] = computeWienerKernel(...
                sumDFFrs, yawAngSpdRS, sampRate, winLen, cutFreq, tauFreq);
            kernels(3).rYawSpd.predResp = computeLinearPrediction(...
                kernels(3).rYawSpd.kernel, sumDFFrs);
            [kernels(3).rYawSpd.varExpl, kernels(3).rYawSpd.fitObj] = ...
                computeVarianceExplained(yawAngSpdRS, ...
                kernels(3).rYawSpd.predResp, 'poly1'); 
            

            % diff kernels, predictions, variance; 4th index
            [kernels(4).fFwdVel.kernel, ~] = computeWienerKernel(...
                fwdVelRS, diffDFFrs, sampRate, winLen, cutFreq, tauFreq);
            kernels(4).fFwdVel.predResp = computeLinearPrediction(...
                kernels(4).fFwdVel.kernel, fwdVelRS);
            [kernels(4).fFwdVel.varExpl, kernels(4).fFwdVel.fitObj] = ...
                computeVarianceExplained(diffDFFrs, ...
                kernels(4).fFwdVel.predResp, 'poly1');               
            
            [kernels(4).rFwdVel.kernel, ~] = computeWienerKernel(...
                diffDFFrs, fwdVelRS, sampRate, winLen, cutFreq, tauFreq); 
            kernels(4).rFwdVel.predResp = computeLinearPrediction(...
                kernels(4).rFwdVel.kernel, diffDFFrs);
            [kernels(4).rFwdVel.varExpl, kernels(4).rFwdVel.fitObj] = ...
                computeVarianceExplained(fwdVelRS, ...
                kernels(4).rFwdVel.predResp, 'poly1');             
            
            [kernels(4).fYawVel.kernel, ~] = computeWienerKernel(...
                yawAngVelRS, diffDFFrs, sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(4).fYawVel.predResp = computeLinearPrediction(...
                kernels(4).fYawVel.kernel, yawAngVelRS);
            [kernels(4).fYawVel.varExpl, kernels(4).fYawVel.fitObj] = ...
                computeVarianceExplained(diffDFFrs, ...
                kernels(4).fYawVel.predResp, 'poly1');             
            
            [kernels(4).rYawVel.kernel, ~] = computeWienerKernel(...
                diffDFFrs, yawAngVelRS, sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(4).rYawVel.predResp = computeLinearPrediction(...
                kernels(4).rYawVel.kernel, diffDFFrs);
            [kernels(4).rYawVel.varExpl, kernels(4).rYawVel.fitObj] = ...
                computeVarianceExplained(yawAngVelRS, ...
                kernels(4).rYawVel.predResp, 'poly1');   
            
            [kernels(4).fYawSpd.kernel, ~] = computeWienerKernel(...
                yawAngSpdRS, diffDFFrs, sampRate, winLen, cutFreq, ...
                tauFreq);
            kernels(4).fYawSpd.predResp = computeLinearPrediction(...
                kernels(4).fYawSpd.kernel, yawAngSpdRS);
            [kernels(4).fYawSpd.varExpl, kernels(4).fYawSpd.fitObj] = ...
                computeVarianceExplained(diffDFFrs, ...
                kernels(4).fYawSpd.predResp, 'poly1');              
            
            [kernels(4).rYawSpd.kernel, ~] = computeWienerKernel(...
                diffDFFrs, yawAngSpdRS, sampRate, winLen, cutFreq, ...
                tauFreq);        
            kernels(4).rYawSpd.predResp = computeLinearPrediction(...
                kernels(4).rYawSpd.kernel, diffDFFrs);
            [kernels(4).rYawSpd.varExpl, kernels(4).rYawSpd.fitObj] = ...
                computeVarianceExplained(yawAngSpdRS, ...
                kernels(4).rYawSpd.predResp, 'poly1');   
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
        
        % save kernel, resampled inputs in struct
        kernelResampInputs.t = newT;
        kernelResampInputs.dFF(1:size(dFFsRS,1),:) = dFFsRS;
        if (numROIs == 2)
            kernelResampInputs.dFF(3,:) = sumDFFrs; % sum as 3rd
            kernelResampInputs.dFF(4,:) = diffDFFrs; % diff as 4th
        end
        kernelResampInputs.fwdVel = fwdVelRS;
        kernelResampInputs.yawAngVel = yawAngVelRS;
        kernelResampInputs.yawAngSpd = yawAngSpdRS;

        
        % save kernels into pData file
%         save('pData.mat', 'kernels', 'kernelParams', ...
%             'kernelResampInputs', '-append');
%         disp('Saved!');
    
        % temporarily, to overwrite previous kernels saved into pData; when
        %  done, replace with above save command
        save('pData.mat', 'kernels', 'kernelParams', ...
            'kernelResampInputs', 'dFFs', 'frameStartTimes', ...
            'numChannels', '-v7.3');
        disp('Saved!');
        
        cd(curDir);
    else
        disp(['Selected trial folder does not contain pDat.mat file' ...
            ' and/or fictracDat.mat']);
        cd(curDir);
        return;
    end
end
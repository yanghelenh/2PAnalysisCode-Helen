% computeSmoFictrac2P.m
%
% Function that takes in one or more pData files and computes a very
%  smoothed version of the FicTrac velocity parameters. To be used for 
%  extracting turning bouts, not to be used as faithful representation of
%  actual velocities.
% Select pData files through GUI. These must have fictrac struct
% Saves smoothed data back into same pData file
% This function is adaptation of
%  LegTrackAnalysis-Helen/computeSmoFictrac.m, except it runs on 2P pData
%
% INPUTS:
%   pDataPath - full path to folder containing pData files
%   sigmaPos - sigma for smoothing position data using smoothdata(), in
%       fictrac samples
%   sigmaVel - sigma for smoothing velocity data using smoothdata(), in
%       fictrac samples
% 
% OUTPUTS:
%   none, but saves fictracSmo output struct back into same pData file
%
% CREATED: 8/25/23
%
% UPDATED: 
%   8/25/23 - HHY
%
function computeSmoFictrac2P(pDataPath, sigmaPos, sigmaVel)

    % default values of sigmaPos and sigmaVel
    if (isempty(sigmaPos))
        sigmaPos = 250;
    end
    if (isempty(sigmaVel))
        sigmaVel = 125;
    end

    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has fictrac, fictracProc, fictracParams
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'fictrac')))
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'fictrac');

        % smooth position
        angPosUn = unwrap(fictrac.yawAngPosWrap);
        angPosUnSmo = smoothdata(angPosUn, 'gaussian',sigmaPos);
        fwdPosSmo = smoothdata(fictrac.fwdCumPos, 'gaussian', sigmaPos);
        slidePosSmo = smoothdata(fictrac.slideCumPos, 'gaussian', sigmaPos);

        % compute velocity
        angVel = gradient(angPosUnSmo, fictrac.t);
        fwdVel = gradient(fwdPosSmo, fictrac.t);
        slideVel = gradient(slidePosSmo, fictrac.t);

        % smooth velocity
        angVelSmo = smoothdata(angVel, 'gaussian',sigmaVel);
        fwdVelSmo = smoothdata(fwdVel, 'gaussian', sigmaVel);
        slideVelSmo = smoothdata(slideVel, 'gaussian', sigmaVel);

        % get bias for yaw and slide velocities - slope of line fit to
        % position over the whole trial
        pAng = polyfit(fictrac.t, angPosUnSmo,1);
        angVelBias = pAng(1);
        pSlide = polyfit(fictrac.t, slidePosSmo, 1);
        slideVelBias = pSlide(1);

        
        % save to output struct
        fictracSmo.yawAngVel = angVelSmo;
        fictracSmo.fwdVel = fwdVelSmo;
        fictracSmo.slideVel = slideVelSmo;

        % add time info to fictracSmo output struct
        fictracSmo.t = fictrac.t;

        % copy over dropInd from fictrac
        fictracSmo.dropInd = fictrac.dropInd;

        % copy over move/not move info from fictrac
        fictracSmo.moveLog = fictrac.moveLog;
        fictracSmo.moveStartInd = fictrac.moveStartInd;
        fictracSmo.moveEndInd = fictrac.moveEndInd;
        fictracSmo.moveStartTimes = fictrac.moveStartTimes;
        fictracSmo.moveEndTimes = fictrac.moveEndTimes;

        % add bias info to fictracSmo output struct
        fictracSmo.angVelBias = angVelBias;
        fictracSmo.slideVelBias = slideVelBias;

        % add sigma info to fictracSmo output struct
        fictracSmo.sigmaPos = sigmaPos;
        fictracSmo.sigmaVel = sigmaVel;

        % append this fictracSmo output struct to same pData file
        save(pDataFullPath, 'fictracSmo', '-append');

        fprintf('Saved fictracSmo for %s!\n', pDataName);

    end
end
% moveCondPairData.m
%
% Function that takes in metaDat struct of selected trials and returns img
%  and fictrac data from pData conditioned on whether the fly was moving or
%  not. Has option to exclude data within specified time window of
%  transition b/w moving and not moving. Interpolates data to specified
%  frame rate. Any data that is excluded (i.e. wrong moving state or too
%  close to transition) is replaced by NaN. Time is not explicitly returned
%  by is implicit in paired data, as no data points are deleted, only
%  replaced by NaNs.
%
% This analyzed data can be used for scatterplots, computing kernels after
%  conditioning on movement, and heat maps.
%
% INPUTS:
%   selMetaDat - metaDat struct with only trials that will be incorporated
%       in this analysis
%   pdPath - full path to pData directory
%   sampRate - sample rate to return conditioned data at
%   moveStartWin - window as [t[before] t[after]] time of transition from
%       not moving to moving that will be excluded from returned paired
%       data, t in seconds
%   moveEndWin - window as [t[before] t[after]] time of transition from
%       moving to not moving that will be excluded from returned paired
%       data, t in seconds
%
% OUTPUTS:
%   condPairData - struct array, 1 element per fly, with fields:
%       img - struct with fields of dF/F for left, right, sum, diff (only
%           those that exist for that fly
%       fictrac - struct with fields of behavioral variables: fwdVel,
%           yawVel, slideVel, yawSpd, totSpd
%       flyID - flyID for which fly contributed that data
%   condPairParams - struct of input parameters
%       sampRate
%       moveStartWin
%       moveEndWin
%   exptNames - names of trials that went into the data
%
% CREATED: 9/12/19 - HHY
% 
% UPDATED:
%   9/12/19 - HHY
%

function [condPairData, condPairParams, exptNames] = moveCondPairData(...
    selMetaDat, pdPath, sampRate, moveStartWin, moveEndWin)

    % fictrac behavioral variables to return
    behVars = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', 'totAngSpd'};

    curDir = pwd;
    cd(pdPath);
    
    % names of all pData files
    pDataFiles = dir(pdPath);
    pDataFileNames = extractfield(pDataFiles, 'name');
    
    % experiment names
    exptNames = selMetaDat.ExperimentName;
    
    % preallocate
    flyIDs = unique(selMetaDat.FlyID);
    numFlies = length(flyIDs);
    
    flyStrct.img.left = [];
    flyStrct.img.right = [];
    flyStrct.img.sum = [];
    flyStrct.img.diff = [];
    
    for k = 1:length(behVars)
        flyStrct.fictrac.move.(behVars{k}) = [];
        flyStrct.fictrac.notMove.(behVars{k}) = [];
    end
    
    flyStrct.flyID = 0;
    flyStrct.moveLog = [];
    flyStrct.notMoveLog = [];
    
    condPairData = repmat(flyStrct, 1, numFlies);
    
    % loop through all flies
    for i = 1:numFlies
        curFlyID = flyIDs(i);
        
        condPairData(i).flyID = curFlyID; % identify fly
        
        % which trials in selMetaDat belong to this fly
        whichInd = find(selMetaDat.FlyID == curFlyID);
        
        % number of trials that belong to this fly
        numTrials = length(whichInd);
        
        % loop through all trials for this fly
        for j = 1:numTrials
            % load relevant pData for this file
            curExptName = selMetaDat.ExperimentName{whichInd(j)};
            pDataName = [curExptName '_pData.mat'];
            
            % check that pData file exists
            hasPDat = sum(strcmp(pDataFileNames, pDataName));
            if (hasPDat)
                % load pData
                load(pDataName, 'fictrac', 'img');
                
                % fictrac's intersample interval, in seconds
                ftIFI = median(diff(fictrac.t));
                
                % number of samples to exclude around move/not move
                % transitions
                moveStartWinSamp = round(moveStartWin ./ ftIFI);
                moveEndWinSamp = round(moveEndWin ./ ftIFI);
                
                % convert edges to logical of excluded time points; 1 for
                %  is too close to transition and should be excluded
                isEdgePoint = false(size(fictrac.moveLog)); % preallocate
                
                % loop through all move start indicies, convert isEdgePoint
                %  logical to 1 at edges
                for k = 1:length(fictrac.moveStartInd)
                    % start of time to exclude (before the transition)
                    beforeInd = fictrac.moveStartInd(k) - ...
                        moveStartWinSamp(1);
                    % account for start of trial
                    if (beforeInd < 1)
                        beforeInd = 1;
                    end
                    % end of time to exclude (after the transition)
                    afterInd = fictrac.moveStartInd(k) + ...
                        moveStartWinSamp(2);
                    % account for end of trial
                    if (afterInd > length(isEdgePoint))
                        afterInd = length(isEdgePoint);
                    end
                    % mark time points to exclude in isEdgePoint
                    isEdgePoint(beforeInd:afterInd) = 1;
                end
                
                % loop through all move end indicies, convert isEdgePoint
                %  logical to 1 at those edges
                for k = 1:length(fictrac.moveEndInd)
                    % start of time to exclude (before the transition)
                    beforeInd = fictrac.moveEndInd(k) - ...
                        moveEndWinSamp(1);
                    % account for start of trial
                    if (beforeInd < 1)
                        beforeInd = 1;
                    end
                    % end of time to exclude (after the transition)
                    afterInd = fictrac.moveEndInd(k) + ...
                        moveEndWinSamp(2);
                    % account for end of trial
                    if (afterInd > length(isEdgePoint))
                        afterInd = length(isEdgePoint);
                    end
                    % mark time points to exclude in isEdgePoint
                    isEdgePoint(beforeInd:afterInd) = 1;
                end
                
                % dropInd as logical
                dropIndLog = false(size(isEdgePoint));
                dropIndLog(fictrac.dropInd) = 1;
                
                % logical 1 for all points to exclude, for moving
                moveNaNsLog = dropIndLog | (~fictrac.moveLog) | ...
                    isEdgePoint;
                % logical 1 for all points to exclude, for not moving
                notMoveNaNsLog = dropIndLog | fictrac.moveLog | ...
                    isEdgePoint;
                    
                % loop through all behavioral variables, generate fictrac
                %  variables with NaNs at appropriate indicies, for moving
                %  and not moving; also removes dropInd                
                for k = 1:length(behVars)
                    tempMove = fictrac.(behVars{k});
                    tempMove(moveNaNsLog) = nan;
                    ftMove.(behVars{k}) = tempMove;
                    
                    tempNotMove = fictrac.(behVars{k});
                    tempNotMove(notMoveNaNsLog) = nan;
                    ftNotMove.(behVars{k}) = tempNotMove;
                end
                
                % interpolate imaging and fictrac data to same sample rate
                newT = img.t(1):(1/sampRate):img.t(end);
                
                % interpolate imaging data
                dffFN = fieldnames(img.filtDFF);
                for k = 1:length(dffFN)
                    dffRS.(dffFN{k}) = interp1(img.t, ...
                        img.filtDFF.(dffFN{k})', newT);
                end
                
                % interpolate fictrac data
                for k = 1:length(behVars)
                    ftMoveRS.(behVars{k}) = interp1(fictrac.t, ...
                        ftMove.(behVars{k})', newT);
                    ftNotMoveRS.(behVars{k}) = interp1(fictrac.t, ...
                        ftNotMove.(behVars{k})', newT);
                end
                
                % interpolate moveNaNsLog, notMoveNaNsLog
                moveNaNsLogRS = logical(interp1(fictrac.t, ...
                    double(moveNaNsLog)', newT));
                notMoveNaNsLog = logical(interp1(fictrac.t, ...
                    double(notMoveNaNsLog)', newT));
                
                % append img data into growing condPairData struct
                imgFN = fieldnames(condPairData(i).img);
                for k = 1:length(imgFN)
                    % if this trial has this ROI type (left, right, sum,
                    %  diff), as not all of them do
                    if (isfield(dffRS, imgFN{k}))
                        condPairData(i).img.(imgFN{k}) = horzcat(...
                            condPairData(i).img.(imgFN{k}), ...
                            dffRS.(imgFN{k}));
                    else % trial does not have this ROI type, append NaNs
                        condPairData(i).img.(imgFN{k}) = horzcat(...
                            condPairData(i).img.(imgFN{k}), ...
                            NaN(size(newT)));
                    end
                    
                    % if this is not the last trial for this fly, add gap
                    %  of 1 NaN to all data (prevents kernel computation
                    %  from using segments from multiple trials)
                    if (j ~= numTrials)
                        condPairData(i).img.(imgFN{k}) = horzcat(...
                            condPairData(i).img.(imgFN{k}), ...
                            NaN(1,1));
                    end
                end
                % append fictrac data into growing condPairData struct
                ftFN = fieldnames(condPairData(i).fictrac.move);
                for k = 1:length(ftFN)
                    % append moving data
                    condPairData(i).fictrac.move.(ftFN{k}) = horzcat(...
                        condPairData(i).fictrac.move.(ftFN{k}), ...
                        ftMoveRS.(ftFN{k}));
                    % append not moving data
                    condPairData(i).fictrac.notMove.(ftFN{k}) = horzcat(...
                        condPairData(i).fictrac.notMove.(ftFN{k}), ...
                        ftNotMoveRS.(ftFN{k}));
                    
                    % if this is not the last trial for fly, add gap of one
                    %  NaN
                    if (j ~= numTrials)
                        condPairData(i).fictrac.move.(ftFN{k}) = ...
                            horzcat(...
                            condPairData(i).fictrac.move.(ftFN{k}), ...
                            NaN(1,1));
                        % append not moving data
                        condPairData(i).fictrac.notMove.(ftFN{k}) = ...
                            horzcat(...
                            condPairData(i).fictrac.notMove.(ftFN{k}), ...
                            NaN(1,1));
                    end
                end 
                % append movement/not movement logicals into growing
                %  condPairData struct
                condPairData(i).moveLog = logical(horzcat(...
                    condPairData(i).moveLog, moveNaNsLogRS));
                condPairData(i).notMoveLog = logical(horzcat(...
                    condPairData(i).notMoveLog, notMoveNaNsLog));
                % if this is not the last trial for fly, add gap of one NaN
                if (j ~= numTrials)
                    condPairData(i).moveLog = logical(horzcat(...
                        condPairData(i).moveLog, true(1,1)));
                    condPairData(i).notMoveLog = logical(horzcat(...
                        condPairData(i).notMoveLog, true(1,1)));
                end
                
            else % if pData file doesn't exist
                fprintf('%s does not exist \n', pDataName);
            end
        end
        
    end
    
    % parameters struct
    condPairParams.sampRate = sampRate;
    condPairParams.moveStartWin = moveStartWin;
    condPairParams.moveEndWin = moveEndWin;
    
    cd(curDir)
end
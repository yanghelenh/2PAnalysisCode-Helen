% preprocessUserDaq.m
%
% Function to take data saved in userDaqDat.mat, the output of run2PExpt
%  for any experiment type, and convert it to a more useable form for later
%  analyses
%
% INPUTS:
%   exptCond - experimental condition string
%   flyData - information about this fly
%   inputParams - experiment-specific settings
%   rawData - data collected on DAQ during experiment
%   rawOutput - signal output on DAQ during experiment
%   settings - static settings about user DAQ
%   sprdshtPath - full path (path/name) of metadata spreadsheet
%   exptName - name of experiment date_fly#_fov#_trial#
%
% OUTPUTS:
%   daqData - data collected on DAQ, with fields labeled
%   daqOutput - signal output on DAQ during experiment, with fields labeled
%   daqTime - timing vector for daqData and daqOutput
%   settings - same as input, passed through
%   
% Also, as side effect, updates metadata spreadsheet with information from
%  flyData and inputParams
%
% CREATED: 10/25/18
% UPDATED: 10/25/18 - HHY
%

function [daqData, daqOutput, daqTime, settings] = preprocessUserDaq(...
    exptCond, flyData, inputParams, rawData, rawOutput, settings,...
    sprdshtPath, exptName)

    % convert rawData array into daqData struct, with named fields
    colNum = 1; % counter of columns in rawData array
    % analog input first
    for i = 1:length(inputParams.aInCh)
        % field names from inputParams
        daqData.(inputParams.aInCh{i}) = rawData(:, colNum);
        colNum = colNum + 1;
    end
    % digital input second
    for i = 1:length(inputParams.dInCh)
        daqData.(inputParams.dInCh{i}) = rawData(:, colNum);
        colNum = colNum + 1;
    end
    
    % convert rawOutput array into daqOutput struct, with name fields
    colNum = 1; % counter of columns in rawOutput array
    % analog input first
    for i = 1:length(inputParams.aOutCh)
        daqOutput.(inputParams.aOutCh{i}) = rawOutput(:, colNum);
        colNum = colNum + 1;
    end
    % digital output second
    for i = 1:length(inputParams.dOutCh)
        daqOutput.(inputParams.dOutCh{i}) = rawOutput(:, colNum);
        colNum = colNum + 1;
    end
    
    % extract timing vector for daqData, daqOutput
    numSamp = size(rawData,1);
    sampRate = settings.bob.sampRate;
    daqTime = (0:(numSamp-1))/sampRate;

    % update metadata spreadsheet
    % read in first column of current metadata spreadsheet, to get number
    %  of rows
    opts = detectImportOptions(sprdshtPath);
    opts.SelectedVariableNames = opts.VariableNames{1};
    sprdshtCol1 = readtable(sprdshtPath, opts);
    % number of rows in spreadsheet currently
    numSsRows = height(sprdshtCol1) + 1; % add 1 for heading row
    
    % starting cell to write new row to
    writeCell = sprintf('A%d', numSsRows + 1);
    
    % build cell array to write to excel spreadsheet
    rowArray = {exptName, exptCond, inputParams.startTimeStamp, ...
        [], flyData.genotype, flyData.manipulation, flyData.prepType, ...
        flyData.age, flyData.ageUnits, inputParams.exptDuration, ...
        inputParams.temperature, [], flyData.dissectionNotes, []};
    
    % write to excel spreadsheet
    xlswrite(sprdshtPath, rowArray, 1, writeCell);
       
end
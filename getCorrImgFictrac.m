% getCorrImgFictrac.m
%
% Function that returns correlation between specified imaging and Fictrac
%  parameters.
% Operates on moveCondData
%
% CREATED: 10/12/23 - HHY
%
% UPDATED:
%   10/12/23 - HHY

function [allFliesCorr, meanCorr, semCorr] = getCorrImgFictrac(...
    whichImg, whichFictrac, datDir)

    % prompt user to select moveCondPairData() files
    [thisFileName, datPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data
    fullFilePath = [datPath filesep thisFileName];
    load(fullFilePath, 'condPairData');

    % number of flies
    numFlies = length(condPairData);

        % preallocate
    allFliesCorr = zeros(numFlies, 1);


    for i = 1:numFlies

        thisFlyImg = condPairData(i).img.(whichImg);
        thisFlyFT = condPairData(i).fictrac.move.(whichFictrac);

        [thisFlyCorr,~] = xcorrWGaps(thisFlyImg, thisFlyFT, 1);

        allFliesCorr(i) = thisFlyCorr;

    end

    % report mean and SEM
    meanCorr = mean(allFliesCorr);
    semCorr = std(allFliesCorr) / sqrt(numFlies);

    fprintf('r = %0.2f +/- %0.2f\n', meanCorr,semCorr);


end
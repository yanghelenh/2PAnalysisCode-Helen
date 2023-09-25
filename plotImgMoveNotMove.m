% plotImgMoveNotMove.m
%
% Function to plot the average dF/F summed b/w the right and left flies
%  during times when the fly is moving compared with times when the fly 
%  is not moving. Uses the output of moveCondPairData(), which computes
%  the moving and not moving periods with transitions excluded.
% Note that moveCondPairData() returns moveLog and notMoveLog as true for
%  points to exclude from named condition
% Select moveCondPairData() output file through GUI
% 
% INPUTS:
%   datDir - path to folder containing moveCondPairData() output files
%
% OUTPUTS:
%
% CREATED: 8/30/23 - HHY
%
% UPDATED:
%   8/30/23 - HHY
%
function p = plotImgMoveNotMove(datDir)
    % prompt user to select moveCondPairData() files
    [thisFileName, datPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data
    fullFilePath = [datPath filesep thisFileName];
    load(fullFilePath, 'condPairData');

    % number of flies
    numFlies = length(condPairData);

    % preallocate (2 for move and not move)
    allFliesAvg = zeros(numFlies, 2);
    allFliesErr = zeros(numFlies, 2);

    for i = 1:numFlies
        % moving, summed dF/F
        thisFlyMoveDat = condPairData(i).img.sum(~condPairData(i).moveLog);
        % not moving, summed dF/F
        thisFlyNotMoveDat = condPairData(i).img.sum(~condPairData(i).notMoveLog);
        
        % get mean for moving
        allFliesAvg(i,1) = mean(thisFlyMoveDat);
        % get SEM for moving
        allFliesErr(i,1) = std(thisFlyMoveDat) / sqrt(length(thisFlyMoveDat));

        % get mean for not moving
        allFliesAvg(i,2) = mean(thisFlyNotMoveDat);
        % get SEM for moving
        allFliesErr(i,2) = std(thisFlyNotMoveDat) / ...
            sqrt(length(thisFlyNotMoveDat));
    end

    % get mean and SEM across flies
    for i = 1:2
        thisCol = allFliesAvg(:,i);
        thisCol(isnan(thisCol)) = [];
        meanAllFlies(i) = mean(thisCol);
        semAllFlies(i) = std(thisCol) / sqrt(length(thisCol));
    end

    % x vector for plotting move and not move 
    xVec = 1:2;
    xVec = xVec';

    figure;
    c = colormap('lines');

    hold on;

    % plot for each individual fly
    plot(xVec, allFliesAvg, ...
        'Marker', '.','LineWidth',0.5, 'Color', c(1,:));

    % plot mean across flies
    errorbar(xVec, meanAllFlies, semAllFlies, ...
        '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    % label x-axis
    xticks(xVec);
    xticklabels({'Moving', 'Not Moving'});

    % label y-axis
    ylabel('Sum dF/F');

    numFlies

    % get stats
    diffAllFlies = allFliesAvg(:,1) - allFliesAvg(:,2);

    [~,p] = ttest(diffAllFlies);
    p

end

% plotDFFCond_bouts.m
%
% Function to plot output of saveDFFCond_bouts(), pooled across
%  multiple flies and conditions. 
% User specifies number of conditions, and then selects files (1 per fly) 
%  for each condition through GUI. Plots mean for each fly and mean across 
%  flies for each condition (different colors for different conditions). 
%  Thin lines for individual flies; thick lines for means across flies.
%
% INPUTS:
%   datDir - full path to directory with output files of
%       saveSpikerateCond_bouts()
%   whichDFF - name of which dF/F value to plot (left, right, sum, diff)
%   numCond - number of conditions to plot (must be matched in num time
%       points)
%   condNames - cell array of length numCond, for names of each condition
%   yScale - 2 element vector for y-axis limits
%
% OUTPUTS:
%   none, but produces plot
%
% CREATED: 8/25/23 - HHY
%
% UPDATED:
%   8/25/23 - HHY
%
function plotDFFCond_bouts(datDir, whichDFF, numCond, condNames)
    
    % preallocate 
    % cell array where each element will be numFlies x numTPts matrix of
    %  mean/SEM for each fly
    allFlyMeans = cell(numCond,1);
    allFlySEMs = cell(numCond,1);

    % cell array where each element will be 1 x numTPts vector of mean/SEM
    %  across flies
    allCondMeans = cell(numCond,1);
    allCondSEMs = cell(numCond,1);

    % loop across number of conditions, get data files for each fly
    % compute mean for each fly and mean across flies
    for i = 1:numCond
        [outputFNames, outputPath] = uigetfile('*.mat', ...
            'Select saveDFFCond_bouts files', ...
            datDir, 'MultiSelect', 'on');

        % if only 1 file selected, not cell array; make sure loop still
        %  works 
        % num flies is number of files
        if (iscell(outputFNames))
            numFlies = length(outputFNames);
        else
            numFlies = 1;
        end

        % preallocate
        thisCondMean = [];
        thisCondSEM = [];

        % loop through all flies
        for j = 1:numFlies
            % handle whether it's a cell array or not
            if (iscell(outputFNames))
                outName = outputFNames{j};
            else
                outName = outputFNames;
            end
            
            outputFullPath = [outputPath outName];

            % load data file
            load(outputFullPath, 'allDFF', 't', 'numBouts');

            thisDFFval = allDFF.(whichDFF);

            thisMean = zeros(size(thisDFFval,2),1);
            thisSEM = zeros(size(thisDFFval,2),1);

            thisSubMean = zeros(50,1);

            % compute mean and SEM for this fly
            for k = 1:size(thisDFFval, 1)
                thisRow = thisDFFval(k,:);
                thisRowMean = mean(thisRow(~isnan(thisRow)));
                thisRowSEM = std(thisRow(~isnan(thisRow))) / ...
                    sqrt(length(thisRow(~isnan(thisRow))));

                thisMean(k) = thisRowMean;
                thisSEM(k) = thisRowSEM;

                if (k>=10 && k<=50)
                    thisSubMean(k-9) = thisRowMean;
                end
            end

%             thisMean = thisMean - mean(thisSubMean);


%             thisMean = mean(thisDFFval, 2); 
%             thisSEM = std(thisDFFval, [], 2) / sqrt(numBouts);

            % save this mean and SEM
            if isrow(thisMean)
                thisMean = thisMean';
            end
            if isrow(thisSEM)
                thisSEM = thisSEM';
            end
            thisCondMean = [thisCondMean, thisMean];
            thisCondSEM = [thisCondSEM, thisSEM];
        end

        % get mean and SEM across flies for this condition
        allCondMeans{i} = mean(thisCondMean,2);
        allCondSEMs{i} = std(thisCondMean, [], 2) / ...
            sqrt(size(thisCondMean,2));

        % save mean and SEM for each fly
        allFlyMeans{i} = thisCondMean;
        allFlySEMs{i} = thisCondSEM;
    end

    % plotting
    figure;
    c = colormap('lines');
    hold on;

    legendInd = [];
    % loop through conditions
    for i = 1:numCond
        % plot individual flies
        plot(t, allFlyMeans{i}, 'LineWidth',0.5, 'Color', c(i,:)');

        % plot mean across flies
        legendInd(i) = plot(t, allCondMeans{i}, ...
            'LineWidth',2, 'Color', c(i,:));
    end

    % plot vertical line at t = 0
    yScale = ylim;

    line([0,0],yScale,'Color','k');

    xlabel('Time relative to yaw peak (s)');
    ylabel('dF/F');
    legend(legendInd,condNames);
end
% plotAvgImgVFictrac.m
%
% Function to plot output of moveCondPairData() as dF/f vs. FicTrac param.
%  Plots average across flies and/or individual flies. Plots data only from
%  moving time points
%
% INPUTS:
%   whichImg - which imaging data to plot ('diff', 'sum', 'right', 'left')
%   whichFictrac - which FicTrac parameter to plot (field of fictrac)
%   datDir - path to folder containing moveCondPairData() output files
%   avg - 'mean' or 'median' for which type of average to plot
%   indivFlies - plot individual flies or not boolean
%   numBins - number of bins in x (FicTrac param)
%   xRange - range of x values to plot and bin across. As 2 element vector 
%       for start and end
%   minBinCount - minimum number of points per fly per bin
%   norm0Vel - boolean for whether to normalize each fly to the value of
%       the bin containing vel = 0
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/29/23 - HHY
%
% UPDATED:
%   8/29/23 - HHY
%
function plotAvgImgVFictrac(whichImg, whichFictrac, datDir, avg, ...
    indivFlies, numBins, xRange, minBinCount, norm0Vel)

    % prompt user to select moveCondPairData() files
    [thisFileName, datPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data
    fullFilePath = [datPath filesep thisFileName];
    load(fullFilePath, 'condPairData');

    % number of flies
    numFlies = length(condPairData);

    % get boundaries of bins
    binSize = (xRange(2) - xRange(1)) / numBins;
    binEdges = xRange(1):binSize:xRange(2);
    binStarts = binEdges(1:(end-1));
    binEnds = binEdges(2:end);
    binMids = (binStarts + binEnds)/2;

    % get index of bin that contains value = 0
    zeroBinInd = find(binStarts <=0 & binEnds > 0);

    % preallocate
    allFliesAvg = zeros(numFlies, length(binMids));
    allFliesErr = zeros(numFlies, length(binMids));


    for i = 1:numFlies

        thisFlyImg = condPairData(i).img.(whichImg);
        thisFlyFT = condPairData(i).fictrac.move.(whichFictrac);

        % loop through all bins, find mean and SEM
        for j = 1:numBins
            % get logical for which samples fall into this bin
            thisBinLog = (thisFlyFT >= binStarts(j)) & ...
                (thisFlyFT < binEnds(j));

            % number of points in this bin
            numPts = sum(thisBinLog);

            if (numPts < minBinCount)
                allFliesAvg(i,j) = nan;
                allFliesErr(i,j) = nan;
            else
                % get img values for these samples
                thisImg = thisFlyImg(thisBinLog);
    
                if (strcmpi(avg,'mean'))
                    % get mean and SEM for this bin
                    allFliesAvg(i,j) = mean(thisImg);
                    allFliesErr(i,j) = std(thisImg) / sqrt(length(thisImg));
                elseif (strcmpi(avg,'median'))
                    % get median and MAD for this bin
                    allFliesAvg(i,j) = median(thisImg);
                    allFliesErr(i,j) = mad(thisImg,1);
                end
            end
        end

        % if we're expressing the curve as relative to the value at vel = 0
        if norm0Vel
            allFliesAvg(i,:) = allFliesAvg(i,:) - allFliesAvg(i,zeroBinInd);
        end
    end



    % generate figure
    figure;

    hold on;

    totAvg = zeros(1,length(binMids));
    totErr = zeros(1, length(binMids));
    % get average and error across flies
    if (strcmpi(avg,'mean'))

        for i = 1:size(allFliesAvg,2)
            thisBin = allFliesAvg(~isnan(allFliesAvg(:,i)),i);
            totAvg(i) = mean(thisBin);
            totErr(i) = std(thisBin) / sqrt(length(thisBin));
        end
    elseif (strcmpi(avg,'median'))
        for i = 1:size(allFliesAvg,2)
            thisBin = allFliesAvg(~isnan(allFliesAvg(:,i)),i);
            totAvg(i) = median(thisBin);
            totErr(i) = mad(thisBin);
        end
    end

    % plot average across flies
%     plot_err_patch_v2(binMids, totAvg, totErr,...
%         [0 0.4470 0.7410],[0.3010 0.7450 0.9330]);

    plot_err_patch_v2(binMids, totAvg, totErr,...
        [0 0 0],[0.5 0.5 0.5]);

    % plot individual flies if flagged, black lines
    if (indivFlies)
        plot(binMids,allFliesAvg');
    end

    % ephys param label
    yLblStr = sprintf('%s (dF/F)', whichImg);
    ylabel(yLblStr);

    xlabel(whichFictrac);
end
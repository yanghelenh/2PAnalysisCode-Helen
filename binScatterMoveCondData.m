% binScatterMoveCondData.m
%
% Function to plot movement conditioned values of one y variable against 1
%  x variable, like a scatterplot, but in x-value bins. For each x bin, 1
%  box and whiskers plot
% Similar to scatterMoveCondData().
% Requires output from moveCondPairData()
% Has option to plot data from all flies on 1 plot or 1 plot per fly
% Has option to plot not moving data as well as moving data. Not moving
%  data is not pooled with moving data; instead, it is plotted with
%  separate box and whiskers at same x-values
% Can specify temporal offset between x and y data. I.e. plot yData
%  against xData from some point in past/future. 
%
% INPUTS:
%   condPairData - array of structs, output of moveCondPairData()
%   xDataName - string specifying name of data to put on x-axis
%   yDataName - string specifying name of data to put on y-axis
%   xScale - x axis limits, as 2 element vector
%   numXBins - number of bins for x variable
%   yScale - y axis limits, as 2 element vector; if empty, matlab
%       automatically scales
%   samePlot - binary for whether to put all flies on same plot (1) or to
%       generate a separate figure for each fly (0)
%   plotNotMove - binary for whether to plot not moving data points; (1)
%       for yes, (0) for no; if yes, plot as x (moving points are circles)
%   offset - offset between xData and yData, positive values are xData
%       before yData, negative values are yData before xData; in units of
%       samples
%   degPerMM - if not empty, will convert any fictrac variables using mm to
%       deg, using this conversion factor
%   ttl - plot title
%
% OUTPUTS:
%   f - handle to figure(s)
%
% CREATED: 9/23/19
%
% UPDATED:
%   9/23/19
%   9/26/19 - added tick marks to x-axis
%   9/27/19 - allow acceleration as fictrac variable, deal with units
%

function f = binScatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, numXBins, yScale, samePlot, plotNotMove, offset, degPerMM, ...
    ttl)

    % fictrac behavioral variables
    behVars = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', ...
        'totAngSpd', 'fwdAcc', 'slideAcc', 'yawAngAcc', 'totAngAccMag'};
    behVarsUnits = {'mm/s', 'mm/s', 'deg/s', 'deg/s', 'deg/s', ...
        'mm/s^2', 'mm/s^2', 'deg/s^2', 'deg/s^2'};
    % imaging variables
    imgVars = {'left', 'right', 'sum', 'diff'};
    
    % number of flies
    numFlies = length(condPairData);
    
    % initalize counter(s) for number of flies that provide valid data
    numFliesMove = 0;
    if (plotNotMove)
        numFliesNotMove = 0;
    end
    
    % define x bins
    % (xmax-xmin)/xNumBins
    xBinWidth = (xScale(2) - xScale(1)) / numXBins;
    % start values of x bins
    xBinStarts = xScale(1):xBinWidth:(xScale(2) - xBinWidth);
    % end values of x bins
    xBinEnds = (xScale(1) + xBinWidth):xBinWidth:xScale(2);
    % midpoints of x bins (for plotting)
    xBinMids = (xBinStarts + xBinEnds) / 2;
    zeroBin = interp1(xBinStarts, 1:length(xBinStarts), 0);

    % preallocate
    allYVals = []; % will be x vector of data of boxplot
    allXVals = []; % will be grouping vector of boxplot
    if (plotNotMove) % for not moving data
        allYValsNM = [];
        allXValsNM = [];
    end
    
    for i = 1:numFlies
        % get x data, assumes xDataName is field of img or fictrac
        if (any(strcmpi(imgVars, xDataName)))
            xDat = condPairData(i).img.(xDataName);
            xDat(condPairData(i).moveLog) = nan;
            % not moving data
            xDatNM = condPairData(i).img.(xDataName);
            xDatNM(condPairData(i).notMoveLog) = nan;
            xUnits = 'dF/F';
        else
            xDat = condPairData(i).fictrac.move.(xDataName);
            xDatNM = condPairData(i).fictrac.notMove.(xDataName);
            % find which behavioral variable it is
            xBehVarInd = find(strcmpi(behVars, xDataName));
            % if the x data is in mm and the user desires a conversion
            %  to degrees
            if ~isempty(degPerMM) && ...
                    strfind(behVarsUnits{xBehVarInd}, 'mm')
                xDat = xDat .* degPerMM;
                xDatNM = xDatNM .* degPerMM;
                xUnits = 'deg/s';
            else
                xUnits = behVarsUnits{xBehVarInd};
            end
        end
        % get y data, assumes yDataName is field of img or fictrac
        if (any(strcmpi(imgVars, yDataName)))
            yDat = condPairData(i).img.(yDataName);
            yDat(condPairData(i).moveLog) = nan;
            % not moving data
            yDatNM = condPairData(i).img.(yDataName);
            yDatNM(condPairData(i).notMoveLog) = nan;
            yUnits = 'dF/F';
        else
            yDat = condPairData(i).fictrac.move.(yDataName);
            yDatNM = condPairData(i).fictrac.notMove.(yDataName);
            % find which behavioral variable it is
            yBehVarInd = find(strcmpi(behVars, yDataName));
            % if the y data is in mm and the user desires a conversion
            %  to degrees
            if ~isempty(degPerMM) && ...
                    strfind(behVarsUnits{yBehVarInd}, 'mm')
                yDat = yDat .* degPerMM;
                yDatNM = yDatNM .* degPerMM;
                yUnits = 'deg/s';
            else
                yUnits = behVarsUnits{yBehVarInd};
            end
        end
        
        % introduce offset to data
        xDat = circshift(xDat, offset);
        yDat = circshift(yDat, offset);
        if (plotNotMove)
            xDatNM = circshift(xDatNM, offset);
            yDatNM = circshift(yDatNM, offset);
        end
        % delete appropriate number of elements from xDat, yDat, zDat, so
        %  wrapped elements aren't used
        if (offset < 0) % negative offsets, remove from end
            xDat = xDat(1:(end + offset));
            yDat = yDat(1:(end + offset));
            if (plotNotMove)
                xDatNM = xDatNM(1:(end + offset));
                yDatNM = yDatNM(1:(end + offset));
            end
        elseif (offset > 0) % positive offsets, remove from front
            xDat = xDat((offset + 1):end);
            yDat = yDat((offset + 1):end);
            if (plotNotMove)
                xDatNM = xDatNM((offset + 1):end);
                yDatNM = yDatNM((offset + 1):end);
            end
        end
        
        % determine if this fly provides any valid data
        % logical array that is true when either xDat or yDat is NaN
        nanLog = isnan(xDat) | isnan(yDat);
        allNaN = all(nanLog); % true if all pairs are NaN
        if (~allNaN) % fly has valid data
            numFliesMove = numFliesMove + 1;
        end
        if (plotNotMove)
            nanLogNM = isnan(xDatNM) | isnan(yDatNM);
            allNanNM = all(nanLogNM);
            if (~allNanNM)
                numFliesNotMove = numFliesNotMove + 1;
            % if fly has no valid data for either movement or not
            % movement, move on to next fly
            elseif (allNaN)
                continue;
            end
        % if not moving data is not being displayed and fly as no valid
        %  moving data, move on to next fly
        elseif (allNaN)
            continue;
        end
        
        % preallocate
        yVals = yDat';
        xVals = zeros(size(yVals));
        yValsNM = yDatNM';
        xValsNM = zeros(size(yValsNM));
        
        % loop through all data points, put each one in appropriate bin
        for j = 1:length(xDat)
            % point has all valid values (not NaN)
            if (~isnan(xDat(j)) && ~isnan(yDat(j)))
                % find index of bin, going from both directions of start
                %  and end of bin
                whichXBinStart = find(xDat(j) >= xBinStarts, 1, 'last');
                whichXBinEnd = find(xDat(j) < xBinEnds, 1, 'first');

                % if either is empty, xDat point exceeds x limits
                if (~isempty(whichXBinStart) && ~isempty(whichXBinEnd))
                    whichXBin = xBinStarts(whichXBinStart);
                else
                    whichXBin = [];
                end

                % note which bin value belongs into
                if (~isempty(whichXBin))
                    xVals(j) = whichXBin;
                    xValsNM(j) = nan;
                    yValsNM(j) = nan;
                else
                    xVals(j) = nan;
                    yVals(j) = nan;
                    xValsNM(j) = nan;
                    yValsNM(j) = nan;
                end
            % if plotting not moving data, if it is valid not moving point
            elseif (plotNotMove && ~isnan(xDatNM(j)) && ~isnan(yDatNM(j)))
                % find index of bin, going from both directions of start
                %  and end of bin
                whichXBinStart = find(xDatNM(j) >= xBinStarts, 1, 'last');
                whichXBinEnd = find(xDatNM(j) < xBinEnds, 1, 'first');

                % if either is empty, xDat point exceeds x limits
                if (~isempty(whichXBinStart) && ~isempty(whichXBinEnd))
                    whichXBin = xBinStarts(whichXBinStart);
                else
                    whichXBin = [];
                end
                
                % note which bin y value belongs into
                if (~isempty(whichXBin))
                    xValsNM(j) = whichXBin;
                    xVals(j) = nan;
                    yVals(j) = nan;
                else
                    xValsNM(j) = nan;
                    yValsNM(j) = nan;
                    xVals(j) = nan;
                    yVals(j) = nan;
                end
            % if point is not valid (has NaN)
            else
                xVals(j) = nan;
                yVals(j) = nan;
                xValsNM(j) = nan;
                yValsNM(j) = nan;
            end
        end
        
        % remove all NaNs, merge into allVals vectors
        yVals(isnan(yVals)) = [];
        xVals(isnan(xVals)) = [];
        
        % last entry is to make sure all xBins appear and so are plotted in
        %  correct x-axis location
        allYVals = [allYVals; yVals; NaN(size(xBinStarts'))];
        allXVals = [allXVals; xVals; xBinStarts'];
        
        if (plotNotMove)
            yValsNM(isnan(yValsNM)) = [];
            xValsNM(isnan(xValsNM)) = [];
            allYValsNM = [allYValsNM; yValsNM; NaN(size(xBinStarts'))];
            allXValsNM = [allXValsNM; xValsNM; xBinStarts'];
        end
        
        % plot, for separate plot per fly
        if (~samePlot)
            f(i) = figure;

            if ((plotNotMove) && ~isempty(allYValsNM) && ~isempty(allXVals))
                boxplot(allYValsNM, allXValsNM, 'PlotStyle', 'compact', ...
                    'Colors', 'r');
            end
            hold on;
            boxplot(allYVals, allXVals, 'PlotStyle', 'compact', ...
                'Colors', 'b');
            
            set(gca, 'xtick', 1:2:length(xBinStarts), ...
                'xticklabels', xBinStarts(1:2:(length(xBinStarts))));
            if ~isempty(yScale)
                ylim(yScale);
            end
            
            xlim([0 numXBins+1]);
            
            xTxt = sprintf('%s (%s)', xDataName, xUnits);
            xlabel(xTxt);
            yTxt = sprintf('%s (%s)', yDataName, yUnits);
            ylabel(yTxt);

            ttlTxt = sprintf('%s %s vs. %s fly ID %d', ttl, ...
                yDataName, xDataName, condPairData(i).flyID);
            title(ttlTxt);
            
            line([0 numXBins+1], [0 0], 'Color', 'k'); % x-axis line
            line([zeroBin zeroBin], yScale, 'Color', 'k'); % y-axis line
            
            % reset for next fly
            allYVals = [];
            allXVals = [];
            if (plotNotMove)
                allYValsNM = [];
                allXValsNM = [];
            end
        end
    end
    
    if (samePlot)
        f = figure;

        if ((plotNotMove) && ~isempty(allYValsNM) && ~isempty(allXVals))
            boxplot(allYValsNM, allXValsNM, 'PlotStyle', 'compact', ...
                'Colors', 'r', 'Whisker', 0, 'Symbol','');
        end
        hold on;
        boxplot(allYVals, allXVals, 'PlotStyle', 'compact', ...
            'Colors', 'b', 'Whisker', 0, 'Symbol','');
        
        set(gca, 'xtick', 1:2:length(xBinStarts), ...
            'xticklabels', xBinStarts(1:2:(length(xBinStarts))));
        
        if ~isempty(yScale)
            ylim(yScale);
        end
        xlim([0 numXBins+1]);

        xTxt = sprintf('%s (%s)', xDataName, xUnits);
        xlabel(xTxt);
        yTxt = sprintf('%s (%s)', yDataName, yUnits);
        ylabel(yTxt);
        
        line([0 numXBins+1], [0 0], 'Color', 'k'); % x-axis line
        line([zeroBin zeroBin], yScale, 'Color', 'k'); % y-axis line

        if (plotNotMove)
            ttlTxt = sprintf('%s %s vs. %s n = %d / %d', ttl, ...
                yDataName, xDataName, numFliesMove, numFliesNotMove);
        else
            ttlTxt = sprintf('%s %s vs. %s n = %d', ttl, ...
                yDataName, xDataName, numFliesMove);
        end
        title(ttlTxt);
    end
end

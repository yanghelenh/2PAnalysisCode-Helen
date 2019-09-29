% scatterMoveCondData.m
%
% Function to generate a scatterplot of movement conditioned data, the
%  output of moveCondPairData().
% Specify data on x-axis and data on y-axis by field name. Can run for any 
%  combination of fictrac or imaging fields.
% Plots multiple flies either on same plot or on different plots
% Plots data during moving bouts but has option to plot data during not
%  moving bouts
% Can specify temporal offset between x and y data. I.e. scatterplot yData
%  against xData from some point in past/future. 
%
% INPUTS:
%   condPairData - array of structs, output of moveCondPairData()
%   xDataName - string specifying name of data to put on x-axis
%   yDataName - string specifying name of data to put on y-axis
%   xScale - x axis limits, as 2 element vector; if empty, matlab
%       automatically scales
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
% CREATED: 9/12/19 - HHY
%
% UPDATED: 9/13/19 - HHY
%   9/27/19 - HHY - handles acceleration, deals with units
%

function f = scatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, yScale, samePlot, plotNotMove, offset, degPerMM, ttl)

    % alpha value for marker transparency
    markerTrans = 0.2;
    % marker size
    markerSize = 20;

    % fictrac behavioral variables
    behVars = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', ...
        'totAngSpd', 'fwdAcc', 'slideAcc', 'yawAngAcc', 'totAngAccMag'};
    behVarsUnits = {'mm/s', 'mm/s', 'deg/s', 'deg/s', 'deg/s', ...
        'mm/s^2', 'mm/s^2', 'deg/s^2', 'deg/s^2'};
    % imaging variables
    imgVars = {'left', 'right', 'sum', 'diff'};
    
    % initalize counter(s) for number of flies that provide valid data
    numFliesMove = 0;
    if (plotNotMove)
        numFliesNotMove = 0;
    end
    
    % if all flies are being plotted on the same plot, initialize figure
    if (samePlot)
        f = figure;
        cmap = colormap('lines');
    else
        numPlots = 1; % counter for number of plots
        g = figure;
        cmap = colormap('lines');
        close(g);
    end
        
    % loop through all flies
    for i = 1:length(condPairData)
        isImg = 0; % variable for whether both x and y are img
        % get x data, assumes xDataName is field of img or fictrac
        if (any(strcmpi(imgVars, xDataName)))
            xDat = condPairData(i).img.(xDataName);
            isImg = 1;
            xUnits = 'dF/F';
        else
            xDat = condPairData(i).fictrac.move.(xDataName);
            % find which behavioral variable it is
            xBehVarInd = find(strcmpi(behVars, xDataName));
            % if the x data is in mm and the user desires a conversion
            %  to degrees
            if ~isempty(degPerMM) && ...
                    strfind(behVarsUnits{xBehVarInd}, 'mm')
                xDat = xDat .* degPerMM;
                xUnits = 'deg/s';
            else
                xUnits = behVarsUnits{xBehVarInd};
            end
        end
        % get y data, assumes yDataName is field of img or fictrac
        if (any(strcmpi(imgVars, yDataName)))
            yDat = condPairData(i).img.(yDataName);
            isImg = isImg * 1;
            yUnits = 'dF/F';
        else
            yDat = condPairData(i).fictrac.move.(yDataName);
            isImg = isImg * 0;
            % find which behavioral variable it is
            yBehVarInd = find(strcmpi(behVars, yDataName));
            % if the y data is in mm and the user desires a conversion
            %  to degrees
            if ~isempty(degPerMM) && ...
                    strfind(behVarsUnits{yBehVarInd}, 'mm')
                yDat = yDat .* degPerMM;
                yUnits = 'deg/s';
            else
                yUnits = behVarsUnits{yBehVarInd};
            end
        end

        % condition x and y data on movement if both are img (haven't
        %  been conditioned on movement yet)
        if (isImg)
            xDat(condPairData(i).moveLog) = nan;
            yDat(condPairData(i).moveLog) = nan;
        end

        % get not move data if it will be plotted
        if (plotNotMove)
            % xDat
            if (sum(strcmpi(imgVars, xDataName)))
                xDatNM = condPairData(i).img.(xDataName);
            else
                xDatNM = condPairData(i).fictrac.notMove.(xDataName);
            end
            % yDat
            if (sum(strcmpi(imgVars, yDataName)))
                yDatNM = condPairData(i).img.(yDataName);
            else
                yDat = condPairData(i).fictrac.notMove.(yDataName);
            end
            % for img data
            if (isImg)
                xDatNM(condPairData(i).notMoveLog) = nan;
                yDatNM(condPairData(i).notMoveLog) = nan;
            end
        end

        % introduce offset to data
        xDat = circshift(xDat, offset);
        if (plotNotMove)
            xDatNM = circshift(xDatNM, offset);
        end
        % delete appropriate number of elements from xDat and yDat, so
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
        
        % scatterplot data
        if (samePlot)
            scatter(xDat(~nanLog), yDat(~nanLog), 'o', 'filled', ...
                'MarkerFaceAlpha', markerTrans, ...
                'MarkerEdgeAlpha', 0, ...
                'SizeData', markerSize, ....
                'MarkerFaceColor', cmap(i, :),...
                'MarkerEdgeColor', cmap(i, :));
            hold on;
            % plot not moving data
            if (plotNotMove && ~allNanNM)
                scatter(xDatNM(~nanLogNM), yDatNM(~nanLogNM), 'x', ...
                    'MarkerFaceAlpha', markerTrans, ...
                    'MarkerEdgeAlpha', markerTrans, ...
                    'SizeData', markerSize, ....
                    'MarkerFaceColor', cmap(i, :),...
                    'MarkerEdgeColor', cmap(i, :));
            end
        else % separate plots for each fly
            f(numPlots) = figure;
            scatter(xDat(~nanLog), yDat(~nanLog), 'o', 'filled', ...
                'MarkerFaceAlpha', markerTrans, ...
                'MarkerEdgeAlpha', 0, ...
                'SizeData', markerSize, ....
                'MarkerFaceColor', cmap(i, :),...
                'MarkerEdgeColor', cmap(i, :));
            hold on;
            % plot not moving data
            if (plotNotMove && ~allNanNM)
                scatter(xDatNM(~nanLogNM), yDatNM(~nanLogNM), 'x', ...
                    'MarkerFaceAlpha', markerTrans, ...
                    'MarkerEdgeAlpha', markerTrans, ...
                    'SizeData', markerSize, ....
                    'MarkerFaceColor', cmap(i+2, :),...
                    'MarkerEdgeColor', cmap(i+2, :));
            end
            
            % add scale, labels for each plot
            if (~isempty(xScale))
                xlim(xScale);
            end
            if (~isempty(yScale))
                ylim(yScale);
            end

            xTxt = sprintf('%s (%s)', xDataName, xUnits);
            xlabel(xTxt);
            yTxt = sprintf('%s (%s)', yDataName, yUnits);
            ylabel(yTxt);

            ttlTxt = sprintf('%s %s vs. %s fly ID %d', ttl, ...
                yDataName, xDataName, condPairData(i).flyID);
            title(ttlTxt);
            
            numPlots = numPlots + 1;
        end
    end
    
    % scale, add labels for when there's one plot
    if (samePlot)
        if (~isempty(xScale))
            xlim(xScale);
        end
        if (~isempty(yScale))
            ylim(yScale);
        end
        
        xTxt = sprintf('%s (%s)', xDataName, xUnits);
        xlabel(xTxt);
        yTxt = sprintf('%s (%s)', yDataName, yUnits);
        ylabel(yTxt);
        
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
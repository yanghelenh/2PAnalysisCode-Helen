close all

xDataName = 'fwdVel';
yDataName = 'sum';
xScale = [-5 25];
numXBins = 30;
yScale = [-1 2];
samePlot = 0;
plotNotMove = 1;
offset = 0;
degPerMM = [];
ttl = 'g34';


f = binScatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, numXBins, yScale, samePlot, plotNotMove, offset, degPerMM, ...
    ttl);
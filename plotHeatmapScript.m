close all

xDataName = 'fwdVel';
yDataName = 'sum';
zDataName = 'counts';
xScale = [-5 25 20];
yScale = [-1 1 20];
zScale = [0 1000];
minNumVals = 20;
offset = 0;
samePlot = 0;
plotNotMove = 1;
degPerMM = [];
ttl = 'g34';

[f, heatmapMat, countsMat] = heatmapMoveCondData(condPairData,...
    xDataName, yDataName, zDataName, xScale, yScale, zScale, minNumVals,...
    offset, samePlot, plotNotMove, degPerMM, ttl);

close all

xDataName = 'yawAngVel';
yDataName = 'fwdVel';
zDataName = 'sum';
xScale = [-500 500 20];
yScale = [-5 25 20];
zScale = [-1 1];
minNumVals = 20;
offset = 0;
samePlot = 0;
plotNotMove = 1;
degPerMM = [];
ttl = 'g34';

[f, heatmapMat, countsMat] = heatmapMoveCondData(condPairData,...
    xDataName, yDataName, zDataName, xScale, yScale, zScale, minNumVals,...
    offset, samePlot, plotNotMove, degPerMM, ttl);

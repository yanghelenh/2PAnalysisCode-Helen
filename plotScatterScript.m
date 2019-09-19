close all

xDataName = 'yawAngVel';
yDataName = 'sum';
xScale = [-500 500];
yScale = [-1 2];
samePlot = 0;
plotNotMove = 1;
offset = 0;
degPerMM = [];
ttl = 'g34';



f = scatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, yScale, samePlot, plotNotMove, offset, degPerMM, ttl);


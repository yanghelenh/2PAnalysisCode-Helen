% plotHeatmapScript.m

close all
clearvars

dataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/230911_moveCond';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190924_scatter-heatmap';

cellType = 'g13';

xDataName = 'yawAngVel';
yDataName = 'fwdVel';
zDataName = 'diff';
xScale = [-450 450 30];
yScale = [-5 25 30];
zScale = [-0.25 0.25];
minNumVals = 20;
offset = 0;
samePlot = 0;
plotNotMove = 1;
degPerMM = [];
ttl = cellType;

% load data for given cell type
load([dataPath filesep 'moveCondDat_' cellType '.mat']);

% plot heat map
[f, heatmapMat, countsMat] = heatmapMoveCondData(condPairData,...
    xDataName, yDataName, zDataName, xScale, yScale, zScale, minNumVals,...
    offset, samePlot, plotNotMove, degPerMM, ttl);

%% save plots
% create folder if it doesn't exist
folderName = sprintf('heatmap_%s_%s v %s and %s', cellType, zDataName, ...
    xDataName, yDataName);
if (~isfolder([figPath filesep folderName]))
    mkdir(figPath, folderName);
end

% save plots
for i = 1:length(f)
    if ~samePlot
        figFileName = sprintf('heatmap_%s_%s v %s %s-fly%d.fig', ...
            cellType, zDataName, xDataName, yDataName, i);
    else
        figFileName = sprintf('heatmap_%s_%s v %s %s-all.fig',...
            cellType, zDataName, xDataName, yDataName);
    end
    fullFigPath = [figPath filesep folderName filesep figFileName];
    saveas(f(i), fullFigPath, 'fig');
end


% plotScatterScript.m

close all
clearvars

dataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190912_moveCond';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190924_scatter-heatmap';

cellType = 'g34';

% parameters
xDataName = 'fwdVel';
yDataName = 'sum';
xScale = [-5 25];
yScale = [-1 2];
samePlot = 0;
plotNotMove = 1;
offset = 0;
degPerMM = [];
ttl = cellType;

% load data for given cell type
load([dataPath filesep 'moveCondDat_' cellType '.mat']);

% generate plot
f = scatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, yScale, samePlot, plotNotMove, offset, degPerMM, ttl);

% save plots
% create folder if it doesn't exist
folderName = sprintf('scatter_%s', cellType);
if (~isfolder([figPath filesep folderName]))
    mkdir(figPath, folderName);
end

% save plots
for i = 1:length(f)
    if ~samePlot
        figFileName = sprintf('scatter_%s-fly%d.fig', cellType, i);
    else
        figFileName = sprintf('scatter_%s-all.fig', cellType);
    end
    fullFigPath = [figPath filesep folderName filesep figFileName];
    saveas(f(i), fullFigPath, 'fig');
end


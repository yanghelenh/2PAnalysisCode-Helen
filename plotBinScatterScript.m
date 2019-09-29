% plotBinScatterScript.m

close all
clearvars

dataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190927_moveCond';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190924_scatter-heatmap';

cellType = 'g13';

% parameters
xDataName = 'totAngSpd';
yDataName = 'sum';
xScale = [0 600];
numXBins = 30;
yScale = [-1 0.5];
samePlot = 1;
plotNotMove = 1;
offset = 0;
degPerMM = [];
ttl = cellType;

% load data for given cell type
load([dataPath filesep 'moveCondDat_' cellType '.mat']);


% generate plots
f = binScatterMoveCondData(condPairData, xDataName, yDataName, ...
    xScale, numXBins, yScale, samePlot, plotNotMove, offset, degPerMM, ...
    ttl);

%% save plots
% create folder if it doesn't exist
folderName = sprintf('binScatter_%s_%s v %s', cellType, yDataName, ...
    xDataName);
if (~isfolder([figPath filesep folderName]))
    mkdir(figPath, folderName);
end

% save plots
for i = 1:length(f)
    if ~samePlot
        figFileName = sprintf('binScatter_%s_%s v %s-fly%d.fig', ...
            cellType, yDataName, xDataName, i);
    else
        figFileName = sprintf('binScatter_%s_%s v %s-all.fig',...
            cellType, yDataName, xDataName);
    end
    fullFigPath = [figPath filesep folderName filesep figFileName];
    saveas(f(i), fullFigPath, 'fig');
end
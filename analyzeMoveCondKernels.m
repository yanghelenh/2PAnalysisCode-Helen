% analyzeMoveCondKernels.m
%
% Script to extract movement conditioned kernels
%
% run extractMoveCondKernels.m on all selected trials, save data, save 
%  figures
%
% NOTE: with current set of movement conditioned data and current kernel
%  parameters (or with winLen = 1.5), the bouts are too short to extract
%  valid kernels
%
% UPDATED: 
%   9/26/19

dataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190912_moveCond';
saveDataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190926_moveCondKernels';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190926_moveCondKernels';

vars = {'Exclude', 'CellType'};

kernelParams.winLen = 3;
kernelParams.cutFreq = 4;
kernelParams.tauFreq = 1;
kernelParams.fwdKernelBW = 5;
kernelParams.revKernelBW = 5;
kernelParams.sampRate = 100;

autoCorrParams.maxLag = 3;

% allCellTypes = {'a01', 'a02', 'b05', 'b06', 'g13', 'g14', 'g15', 'g16', ...
%     'g31', 'g34', 'p05', 'p09', 'p11', 'p12', 'p18', 'p32'};
allCellTypes = {'g13'};

% load y axis scale and labels; variables: yScale, yLabels, yScaleAllDeg,
% yLabelsAllDeg
loadKernelYScaleYLabels();

acYScale = {[-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1],...
    [-0.5 1], [-0.5 1], [-0.5 1]};

%%


% loop through all cell types
for i = 1:length(allCellTypes)
    disp(allCellTypes{i});
    conds = {'~=1', allCellTypes{i}};
    
    % load data
    load([dataPath filesep 'moveCondDat_' allCellTypes{i} '.mat']);
    
    % compute all kernels, this takes a while    
    [kernels, kernelParams, autoCorrParams, autoCorr] = ...
        extractMoveCondKernels(condPairData, kernelParams, autoCorrParams);
    
    % save data
    save([saveDataPath filesep allCellTypes{i} '.mat'], ...
        'kernels', 'autoCorr', 'kernelParams', 'autoCorrParams', ...
        'vars', 'conds', '-v7.3');
    
    % generate kernel figures
    
    kernelSEMFig = plotMeanKernels(kernels, kernelParams, yScaleAllDeg, ...
        allCellTypes{i}, yLabelsAllDeg, degPerMM, 1, 0, 0);
    kernelIndivFig = plotMeanKernels(kernels, kernelParams, yScaleAllDeg, ...
        allCellTypes{i}, yLabelsAllDeg, degPerMM, 0, 1, 0);
    kernelMeanOnlyFig = plotMeanKernels(kernels, kernelParams, yScaleAllDeg, ...
        allCellTypes{i}, yLabelsAllDeg, degPerMM, 0, 0, 0);
    
    % generate autocorrelation figures
    autoCorrSEMFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
        allCellTypes{i}, 1, 0);
    autoCorrIndivFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
        allCellTypes{i}, 0, 1);
    autoCorrMeanFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
        allCellTypes{i}, 0, 0);
    
    
    
    % save figures
    saveas(kernelSEMFig, ...
        [figPath filesep allCellTypes{i} '_kernels_semFig'], 'fig');
    saveas(kernelIndivFig, ...
        [figPath filesep allCellTypes{i} '_kernels_indivFig'], 'fig');
    saveas(kernelMeanOnlyFig, ...
        [figPath filesep allCellTypes{i} '_kernels_meanOnlyFig'], 'fig');
    
    saveas(autoCorrSEMFig, ...
        [figPath filesep allCellTypes{i} '_autoCorr_semFig'], 'fig');
    saveas(autoCorrIndivFig, ...
        [figPath filesep allCellTypes{i} '_autoCorr_indivFig'], 'fig');
    saveas(autoCorrMeanFig, ...
        [figPath filesep allCellTypes{i} '_autoCorr_meanOnlyFig'], 'fig');
    
    % close all figures
%     close all
    
end
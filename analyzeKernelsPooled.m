% analyzeKernelsPooled.m
%
% Script to compute kernels averaged over multiple trials. Saves the
% kernels as well as the corresponding figures
%
% run extractKernels.m on all selected trials, save data, save figures
%
% UPDATED: 
%   8/29/19
%   9/4/19 - autocorrelation computation now working

datPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190904';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190904_kernels-autoCorr';

vars = {'Exclude', 'CellType'};

kernelParams.winLen = 3;
kernelParams.cutFreq = 4;
kernelParams.tauFreq = 1;
kernelParams.fwdKernelBW = 5;
kernelParams.revKernelBW = 5;
kernelParams.sampRate = 100;

autoCorrParams.maxLag = 3;

allCellTypes = {'a01', 'a02', 'b05', 'b06', 'g13', 'g14', 'g15', 'g16', ...
    'g31', 'g34', 'p05', 'p09', 'p11', 'p12', 'p18', 'p32'};

% load y axis scale and labels; variables: yScale, yLabels, yScaleAllDeg,
% yLabelsAllDeg
loadKernelYScaleYLabels();

acYScale = {[-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1],...
    [-0.5 1], [-0.5 1], [-0.5 1]};

%%
% load metadata spreadsheet
metaDat = loadMetadataSpreadsheet();


% loop through all cell types
for i = 1:length(allCellTypes)
    disp(allCellTypes{i});
    conds = {'~=1', allCellTypes{i}};
    
    % select appropriate pData
    [~, selMetaDat] = returnSelectMetaDat(metaDat, vars, conds);

    % compute all kernels, this takes a while
    [kernels, autoCorr, exptNames, kernelParams, autoCorrParams] = ...
        extractKernels(...
        selMetaDat, pDataPath(), kernelParams, autoCorrParams);
    
    % save data
    save([datPath filesep allCellTypes{i} '.mat'], 'selMetaDat', ...
        'kernels', 'autoCorr', 'exptNames', 'kernelParams', ...
        'autoCorrParams', 'vars', 'conds', '-v7.3');
    
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
    close all
    
end
% analyzeKernelsPooled.m
%
% Script to compute kernels averaged over multiple trials. Saves the
% kernels as well as the corresponding figures
%
% run extractKernels.m on all selected trials, save data, save figures
%
% UPDATED: 8/29/19

datPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190828';
figPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/Figures/190828_kernels';

vars = {'Exclude', 'CellType'};

kernelParams.winLen = 3;
kernelParams.cutFreq = 4;
kernelParams.tauFreq = 1;
kernelParams.fwdKernelBW = 5;
kernelParams.revKernelBW = 5;
kernelParams.sampRate = 100;

allCellTypes = {'a01', 'a02', 'b05', 'b06', 'g13', 'g14', 'g15', 'g16', ...
    'g31', 'g34', 'p05', 'p09', 'p11', 'p12', 'p18', 'p32'};

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
    [kernels, ~, exptNames, kernelParams] = extractKernels(...
        selMetaDat, pDataPath(), kernelParams);
    
    % save data
    save([datPath filesep allCellTypes{i} '.mat'], 'selMetaDat', ...
        'kernels', 'exptNames', 'kernelParams', 'vars', 'conds', '-v7.3');
    
    % generate figures
    semFig = plotMeanKernels(kernels, kernelParams, [], ...
        allCellTypes{i}, 1, 0);
    indivFig = plotMeanKernels(kernels, kernelParams, [], ...
        allCellTypes{i}, 0, 1);
    meanOnlyFig = plotMeanKernels(kernels, kernelParams, [], ...
        allCellTypes{i}, 0, 0);
    
    % save figures
    saveas(semFig, [figPath filesep allCellTypes{i} '_semFig'], 'fig');
    saveas(indivFig, [figPath filesep allCellTypes{i} '_indivFig'], 'fig');
    saveas(meanOnlyFig, [figPath filesep allCellTypes{i} '_meanOnlyFig'], ...
        'fig');
    
    % close all figures
    close all
    
end
% saveAllMoveCondData.m
%
% Script to generate and save all movement conditioned, paired imaging and
%  fictrac data for flies that meet selection criteria
%
% Selects flies/trials using metadata spreadsheet and returnSelectMetaDat()
% Runs moveCondPairData() on all selected trials and saves data
%
% UPDATED: 
%   9/13/19
%

datPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen/AnalyzedData/190912_moveCond';

vars = {'Exclude', 'CellType'};

allCellTypes = {'a01', 'a02', 'b05', 'b06', 'g13', 'g14', 'g15', 'g16', ...
    'g31', 'g34', 'p05', 'p09', 'p11', 'p12', 'p18', 'p32'};

% parameters for computing moveStart and moveEnd windows
tPre = 0.25; % how much calcium precedes behavior, in sec
tauOn = 0.07; % half rise time for 10 APs for jGCaMP7s
tauOff = 1.69; % half decay time for 10 APs for jGCaMP7s


% define move start and move end windows using above parameters, assumes
%  that calcium signal is positively correlated with movement
moveStartWin = [tPre (2*tauOn)];
moveEndWin = [tPre (2*tauOff)];

sampRate = 100; % 100 Hz, like in kernel computations

%%
% load metadata spreadsheet
metaDat = loadMetadataSpreadsheet();

% loop through all cell types
for i = 1:length(allCellTypes)
    disp(allCellTypes{i});
    conds = {'~=1', allCellTypes{i}};
    
    % select appropriate pData
    [~, selMetaDat] = returnSelectMetaDat(metaDat, vars, conds);

    % compute movement conditioned data
    [condPairData, condPairParams, exptNames] = moveCondPairData(...
        selMetaDat, pDataPath(), sampRate, moveStartWin, moveEndWin);

    % save data
    save([datPath filesep 'moveCondDat_' allCellTypes{i} '.mat'], ...
        'selMetaDat', 'condPairData', 'condPairParams', 'exptNames', ...
        'vars', 'conds', '-v7.3');
    
end
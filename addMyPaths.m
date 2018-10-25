% addMyPaths.m
% 
% script to add relevant 2p data analysis paths
%
% Users: change the following paths to match those on your local computer

%% ScanImageTiffReader code
scanImageTiffReaderPath = ...
    '/Users/hyang/Documents/MATLAB/ScanImageTiffReader';
addpath(scanImageTiffReaderPath);

%% Animal Part Tracker code
APTpath = '/Users/hyang/Documents/MATLAB/APT';
addpath(APTpath);
    
%% Analysis code repository (2PAnalysisCode-Helen)
analysisPath = '/Users/hyang/Documents/2PAnalysisCode-Helen';
addpath(analysisPath);
    
%% Experimental code repository (2PCode-Helen)
expPath = '/Users/hyang/Documents/2PCode-Helen';
addpath(expPath);
    
%% Folder containing your metadata spreadsheet 
metadataPath = '/Users/hyang/Documents/2PAnalysis-Helen';
addpath(metadataPath);
    
%% Folder containing processed data (*_pdata.mat files)
% pDataPath = '/Users/hyang/Documents/2PAnalysis-Helen/pData';
% addpath(pDataPath);
% addMyPaths.m
% 
% function to add relevant 2p data analysis paths, including subfolders
%
% Users: change the following paths to match those on your local computer

function addMyPaths()
    %% ScanImageTiffReader code
    scanImageTiffReaderPath = ...
        '/Users/hyang/Documents/MATLAB/ScanImageTiffReader';
    addpath(genpath(scanImageTiffReaderPath));

    %% Animal Part Tracker code
    APTpath = '/Users/hyang/Documents/MATLAB/APT';
    addpath(genpath(APTpath));

    %% Analysis code repository (2PAnalysisCode-Helen)
    analysisPath = '/Users/hyang/Documents/2PAnalysisCode-Helen';
    addpath(genpath(analysisPath));

    %% Experimental code repository (2PCode-Helen)
    expPath = '/Users/hyang/Documents/2PCode-Helen';
    addpath(genpath(expPath));

    %% Folder containing your metadata spreadsheet 
%     metadataPath = '/Users/hyang/Documents/2PAnalysis-Helen';
    metadataPath = '/Users/hyang/Dropbox (HMS)/2PAnalysis-Helen';
    addpath(genpath(metadataPath));
    
    %% Python functions, add to python system path
    a2libPath = '/Users/hyang/Documents/2PAnalysisCode-Helen/a2lib';
    P = py.sys.path;
    if count(P,a2libPath) == 0
        insert(P,int32(0),a2libPath);
    end

    %% Folder containing processed data (*_pdata.mat files)
    addpath(pDataPath());

end
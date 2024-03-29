% preprocessFicTracOnly.m
%
% Temporary variant of preprocess just to update FicTrac data (because
%  previous lowpas filter cutoff wrong).
%
% On imaging data:
% Reads in ScanImage .tif files. Aligns each time series. Saves the aligned
%  series, unaligned series, mean image, and metadata output by ScanImage.
% 
% On FicTrac data:
% Converts raw voltage traces collected on experimental DAQ into position
%  and velocity traces for each of 3 axes (yaw/heading, forward, slide).
%  Also generates fictive path. Saves all of these in .mat file.
% 
% On leg tracking data:
% Extracts frame time for each frame of leg tracking video. Saves this as
%  well as path and name info for legVid video, which is to be fed into
%  APT.
%
% NOTE: assumes folder organization of date folder with one or more fly
% folders with one or more field of view folders with one or more trial
% folders. Assumes trial folder only has 1 f*.tif file, which is ScanImage
% .tif file. ScanImage .tif file name must begin with f (i.e. default 
% starts with 'file'. Assumes trial folder only has 1 *legVid*.mp4 file, 
% which is the leg vid video file.
%
% CREATED: 10/3/18 HHY
% UPDATED: 4/15/19 HHY 
%
function preprocessFicTracOnly()

    disp('Select metadata spreadsheet for experiment');
    [sprdshtName, sprdshtPath] = uigetfile('*.xlsx', ...
        'Select metadata spreadsheet');
    sprdshtFullPath = [sprdshtPath filesep sprdshtName];

    disp('Select a date folder to preprocess.');
    datePath = uigetdir;
    
    % get name of date folder
    fsLoc = strfind(datePath, filesep);
    dateFolder = datePath((fsLoc(end) + 1):end);
    
    curDir = pwd;
    cd(datePath);
    
    
    % get all fly folders in date directory    
    flyFolders = dir([datePath filesep 'fly*']);

    % loop through folders in date folder
    for i = 1:length(flyFolders)
        flyPath = [datePath filesep flyFolders(i).name];
        
        % get all FOV folders in fly folder
        fovFolders = dir([flyPath filesep 'fov*']);
        
        % loop through folders in fly folder
        for j = 1:length(fovFolders)
            fovPath = [flyPath filesep fovFolders(j).name];
            
            % get all trial folders in FOV folder
            trialFolders = dir([fovPath filesep 'trial*']);
            
            % loop through trial folders in FOV folder
            for k = 1:length(trialFolders)
                trialPath = [fovPath filesep trialFolders(k).name];
                exptName = [dateFolder '_' flyFolders(i).name '_' ...
                    fovFolders(j).name '_' trialFolders(k).name];
                cd (trialPath)
                
                % load data from experimental DAQ, metadata
                try
                    load([trialPath filesep 'userDaqDat.mat']);
                catch mExcep
                    % if there is no userDaqDat.mat file
                    if(strcmp(mExcep.identifier, ...
                            'MATLAB:load:couldNotReadFile'))
                        fprintf(['Error: no userDaqDat.mat file found '...
                            'for \n %s \n Skipping processing'], ... 
                            trialPath);
                        continue; % skips processing of this trial
                    % other errors, throw and stop preprocess completely    
                    else 
                        rethrow(ME); 
                    end
                end
                
                % display updates in command line
                fprintf('Preprocessing %s \n', trialPath);
                
                % process metadata from userDaqDat.mat
                [daqData, daqOutput, daqTime, settings] = ...
                    preprocessUserDaq(exptCond, flyData, inputParams, ...
                    rawData, rawOutput, settings, sprdshtFullPath,...
                    exptName);
                
               
                
                % if this experiment has FicTrac data
                if (contains(exptCond, 'Fictrac'))
                    % Process FicTrac data
                    disp('Preprocessing FicTrac data');
                    preprocessFicTrac(daqData, daqTime, ...
                        settings.bob.sampRate);
                end
                
                
                % display updates in command line
                fprintf('Done preprocessing %s \n', trialPath);
                
            end
        end
    end
    
    cd(curDir);
end
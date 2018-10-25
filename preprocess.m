% preprocess.m
%
% First function to run on imaging/FicTrac/leg tracking data. 
%
% On imaging data:
% Reads in ScanImage .tif files. Aligns each time series. Saves the aligned
%  series, unaligned series, mean image, and metadata output by ScanImage.
% 
% On FicTrac data:
% 
% On leg tracking data:
%
% NOTE: assumes folder organization of date folder with one or more fly
% folders with one or more field of view folers with one or more trial
% folders. Assumes trial folder only has 1 .tif file, which is ScanImage
% .tif file
%
% CREATED: 10/3/18 HHY
% UPDATED: 10/25/18 HHY 
%
function preprocess()

    disp('Select metadata spreadsheet for experiment');
    [sprdshtName, sprdshtPath] = uigetfile;
    sprdshtFullPath = [sprdshtName filesep sprdshtPath];

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
                expName = [dateFolder '_' flyFolders(i).name '_' ...
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
                            'for /n %s /n Skipping processing'], ... 
                            trialPath);
                        continue; % skips processing of this trial
                    % other errors, throw and stop preprocess completely    
                    else 
                        rethrow(ME); 
                    end
                end
                
                % process metadata from userDaqDat.mat
                [daqData, daqOutput, daqTime, settings] = ...
                    preprocessUserDaq(exptCond, flyData, inputParams, ...
                    rawData, rawOutput, settings, sprdshtPath, exptName);
                
                
                %% if this experiment has imaging data
                if (contains(exptCond, 'Img'))
                    % name of ScanImage Tiff
                    tifFile = dir([trialPath filesep '*.tif']);

                    % only if imaging data exists
                    if ~isempty(tifFile)
                        % align imaging data, save that and metadata in
                        %  trialPath
                        preprocessImaging(tifFile, daqData, daqTime);
                    else
                        fprintf(['Warning: Imaging data expected, but '
                            'no .tif file found for:/n %s'], trialPath);
                    end
                end
                
                %% if this experiment has FicTrac data
                if (contains(exptCond, 'Fictrac'))
                    % Process FicTrac data
                end
                
                %% if this experiment has leg tracking data
                if (contains(exptCond, 'leg'))
                    % Process leg data
                end
                
            end
        end
    end
    
    cd(curDir);
end
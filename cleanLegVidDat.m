% cleanLegVidDat.m
%
% Quick function to go through all date folders preprocessed on or before
%  2/25/19 to fix error in getting vidName. Was saving hidden files like
%  ._legVid.mp4 as the legVid name instead of legVid.mp4. Fixed in
%  preprocessLegVid.m on 2/25/19
%
% CREATED: 2/25/19 HHY
%


function cleanLegVidDat()    
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
                    load([trialPath filesep 'legVidDat.mat']);
                catch mExcep
                    % if there is no userDaqDat.mat file
                    if(strcmp(mExcep.identifier, ...
                            'MATLAB:load:couldNotReadFile'))
                        fprintf(['Error: no legVidDat.mat file found '...
                            'for \n %s \n Skipping processing'], ... 
                            trialPath);
                        continue; % skips processing of this trial
                    % other errors, throw and stop preprocess completely    
                    else 
                        rethrow(ME); 
                    end
                end
                
                vidDir = dir('legVid*.mp4');

                % check that legVid video exists
                if (~isempty(vidDir))
                    vidName = vidDir.name;
                    vidFolder = vidDir.folder;
                    vidExists = true;
                else
                    disp('Warning: No leg video matching name legVid*.mp4 found');
                    vidName = [];
                    vidFolder = [];
                    vidExists = false;
                end

                save('legVidDat.mat', 'legVidFrameTimes', 'vidName', ...
                    'vidFolder','vidExists', '-v7.3');
               
                
            end
        end
    end
    
%     cd(curDir);
end
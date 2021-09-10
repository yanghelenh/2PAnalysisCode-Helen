% moveLegVideos.m
%
% Function to move all leg videos in input date directories to a 
%  user-selected destination folder (through GUI). Renames leg videos with
%  full experiment name (date_fly#_fov#_trial#_legVid.mp4).
%
% INPUT:
%   dateDirs - cell array of full paths to date directories
%
% OUTPUT:
%   none, but copies all leg videos to user selected folder with experiment
%       name
%
% CREATED: 9/30/19 - HHY
%
% UPDATED: 
%   9/30/19 - HHY
%

function moveLegVideos(dateDirs)

    disp('Select destination folder for videos.');
    vidPath = uigetdir;

    for d = 1:length(dateDirs)

        datePath = dateDirs{d};
        
        % get name of date folder
        fsLoc = strfind(datePath, filesep);
        dateFolder = datePath((fsLoc(end) + 1):end);

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
                    
                    % if video file exists
                    if (isfile('legVid.mp4'))
                        newFilePath = [vidPath filesep exptName '_legVid.mp4'];
                        copyfile('legVid.mp4', newFilePath);
                        
                    end
                end
            end
        end
    end

end
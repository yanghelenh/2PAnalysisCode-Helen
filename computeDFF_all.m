% computeDFF_all.m
%
% Function to compute dF/F for all date folders specified in function input
%  or through gui.
%
% INPUT:
%   dateFolderPaths - [] if gui to select one date folder; otherwise,
%       cell array of full paths to date folders
%
% OUTPUTS:
%   None, but creates pData.mat file to save dF/F, timing, and # of
%     channels info
%
% CREATED: 7/30/19
% UPDATED: 7/30/19
%   8/27/19 - OBSOLETE - works but is no longer in analysis pipeline
%

function computeDFF_all(dateFolderPaths)

    % user selects one date folder if dateFolderPaths isn't a cell array of
    %  strings already
    if ~iscellstr(dateFolderPaths)
    
        disp('Select a date folder to get kernels for.');
        datePath = uigetdir;
        dateFolderPaths = {datePath};
    end

    curDir = pwd;
    
    
    for d = 1:length(dateFolderPaths)
        datePath = dateFolderPaths{d};
        
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
                    cd (trialPath)

                    computeDFF_trial(trialPath);

                end
            end
        end
        fprintf('Done computing dFF for \n %s \n', datePath');
        cd(curDir);
    end 
end

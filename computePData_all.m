% computePData_all.m
%
% Function to compute all processed data (currently, for imaging and
%  FicTrac data) for all trials in specified date folder(s). 
% Saves processed data into a pData.mat file in the same trial folder.
% Preprocessing, ROI selection, and FicTrac dropped frames selection must
%  have been completed
% Calls computePData_trial.m
% Note that if pData.mat file already exists, this function overwrites it.
%
% INPUT:
%   dateFolderPaths - [] if gui to select one date folder; otherwise,
%       cell array of full paths to date folders
%
% OUTPUT:
%   none, but creates pData.mat file for all trials
%
% CREATED: 8/28/19 - HHY
% UPDATED: 8/28/19 - HHY
%

function computePData_all(dateFolderPaths)

    % user selects one date folder if dateFolderPaths isn't a cell array of
    %  strings already
    if ~iscellstr(dateFolderPaths)
    
        disp('Select a date folder to get pData for.');
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

                    % compute pData for trial
                    computePData_trial(trialPath);

                end
            end
        end
        fprintf('Done computing pData for \n %s \n', datePath');
        cd(curDir);
    end 
    
    cd(curDir);
end

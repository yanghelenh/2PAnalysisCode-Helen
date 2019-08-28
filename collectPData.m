% collectPData.m
%
% Function to copy all pData in date folder(s) to a pData folder. Renames
%  pData file to include experiment name in front.
%
% INPUT:
%   dateFolderPaths - [] if gui to select one date folder; otherwise,
%       cell array of full paths to date folders
%   pdPath - full path to pData folder
%
% OUTPUT:
%   none, but copies all pData from date folder(s) to a pData folder
%
% CREATED: 8/28/19 - HHY
% UPDATED: 8/28/19 - HHY
%

function collectPData(dateFolderPaths, pdPath)

    % user selects one date folder if dateFolderPaths isn't a cell array of
    %  strings already
    if ~iscellstr(dateFolderPaths)
    
        disp('Select a date folder to collect pData for.');
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

                    % collect pData for trial
                    % get file info for pData file
                    pDataFile = dir('pData.mat');
                    if ~isempty(pDataFile)
                        if ~(isempty(who('-file', 'pData.mat', 'name')))
                            load('pData.mat', 'name');
                            copyfile(pDataFile.name, ...
                                [pdPath filesep ...
                                name.exptName '_pData.mat']);
                        else
                            fprintf(...
                                '%s \n contains old pData.mat file \n',...
                                trialPath);
                        end
                    else
                    	fprintf(...
                            '%s \n does not contain pData.mat file \n',...
                            trialPath);
                    end

                end
            end
        end
        fprintf('Done collecting pData for \n %s \n', datePath');
        cd(curDir);
    end 
    
    cd(curDir);
    
end
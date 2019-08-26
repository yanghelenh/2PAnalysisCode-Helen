% plotAllKernelsDateDir()
%
% Function to plot all kernels given either cell array of paths to date
% directories or uses user selected single date dir.
%
% INPUT:
%   dateFolderPaths - [] if gui to select one date folder; otherwise,
%       cell array of full paths to date folders
%
% OUTPUT:
%   none but generates plots
%
% CREATED: 5/24/19
% UPDATED: 5/24/19
%

function plotAllKernelsDateDir(dateFolderPaths)

    % user selects one date folder if dateFolderPaths isn't a cell array of
    %  strings already
    if ~iscellstr(dateFolderPaths)
    
        disp('Select a date folder to plot kernels for.');
        datePath = uigetdir;
    
        dateFolderPaths = {datePath};
    end

    curDir = pwd;

    for d = 1:length(dateFolderPaths)
        datePath = dateFolderPaths{d};
        fsLoc = strfind(datePath, filesep);
        dateFolder = datePath((fsLoc(end) + 1):end);
        
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
                    
                    exptName = [dateFolder '\_' flyFolders(i).name '\_' ...
                        fovFolders(j).name '\_' trialFolders(k).name];

                    % check that pData file exists
                    trialPathFiles = dir(trialPath);
                    trialPathFileNames = extractfield(trialPathFiles, ...
                        'name');
                    hasPDat = sum(strcmp(trialPathFileNames, 'pData.mat'));
                    
                    if hasPDat
                        if (isempty(who('-file', 'pData.mat', 'kernels')))
                            fprintf('%s lacks kernels \n', trialPath);
                        else
                            load('pData.mat', 'kernels', ...
                                'kernelParams');
                            plotAllKernels(kernels, kernelParams, ...
                                exptName);
                        end
                    end

                end
            end
        end
        fprintf('Done plotting %s \n', datePath');
        cd(curDir);
    end 
end
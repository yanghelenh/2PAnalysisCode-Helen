% saveAllKernels.m
%
% Function to compute kernels for all selected date directories.
% Calls extractAllKernels()
%
% INPUT:
%   dateFolderPaths - [] if gui to select one date folder; otherwise,
%       cell array of full paths to date folders
%   winLen - length of window, in seconds that averages are computed over
%   cutFreq - cutoff frequency (f_cut) in attenuation applied to frequency
%       domain filter, a la Nagel and Wilson 2011. Set to 0 with tauFreq if
%       no attenuation desired.
%   tauFreq - f_tau in attenuation applied to frequency domain filter, a la
%       Nagel and Wilson 2011. Set to 0 with cutFreq if no attenuation
%       desired.
%   sampRate - sampling rate to convert dF/F and FicTrac data to, and to
%       calculate kernel at
%
% OUTPUT:
%   none - but saves kernels in pData.mat files
%
% CREATED: 5/24/19
% UPDATED: 5/24/19
%

function saveAllKernels(dateFolderPaths, winLen, cutFreq, tauFreq, ...
    sampRate)

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

                    extractAllKernels(trialPath, winLen, cutFreq, ...
                        tauFreq, sampRate, 0);

                end
            end
        end
        fprintf('Done computing kernels for \n %s \n', datePath');
        cd(curDir);
    end 
end

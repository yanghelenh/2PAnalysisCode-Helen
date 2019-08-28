% getExptName.m
%
% Helper function called by computePData_trial.m that turns full trial path
%  into an experiment name struct, with fields for the date folder, fly,
%  fov, and trial, as well as for the full name of the experiment.
%
% INPUT:
%   trialPath - full trial path, as string
%
% OUTPUT:
%   name - struct containing name info for trial, all strings
%       dateName
%       flyName
%       fovName
%       trialName
%       exptName - full experiment name DATE_fly##_fov##_trial##
%
% CREATED: 8/28/19 - HHY
%
% UPDATED: 8/28/19 - HHY
%

function name = getExptName(trialPath)

    % find file separators
    sepInd = strfind(trialPath, filesep);
    
    % trial name is last file separator to end
    trialName = trialPath((sepInd(end)+1):end);
    
    % fov name is between last 2 file separators
    fovName = trialPath((sepInd(end - 1)+1):(sepInd(end)-1));
    
    % fly name is between 3rd to last and 2nd to last file separators
    flyName = trialPath((sepInd(end - 2)+1):(sepInd(end - 1)-1));
    
    % date name is between 4th to last and 3rd to last file separators
    dateName = trialPath((sepInd(end - 3)+1):(sepInd(end - 2)-1));
    
    % construct full experiment name
    exptName = [dateName '_' flyName '_' fovName '_' trialName];
    
    % add all tro name struct
    name.dateName = dateName;
    name.flyName = flyName;
    name.fovName = fovName;
    name.trialName = trialName;
    name.exptName = exptName;
end
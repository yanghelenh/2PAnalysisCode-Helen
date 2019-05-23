% userSelectTrial.m
%
% Function that brings up the user interface for them to select the trial
%  folder to process.
%
% INPUTS:
%   none - but prompts user through folder-selection gui
%
% OUTPUTS:
%   uTrialPath - full path of selected trial folder
%   curDir - present working directory
%
% CREATED: 5/21/19
% UPDATED: 5/21/19
%

function [uTrialPath, curDir] = userSelectTrial()
    % ask user to select trial folder
    disp('Select a trial folder to display.');
    uTrialPath = uigetdir;
    curDir = pwd;
    
    fprintf('%s \n', uTrialPath);
end
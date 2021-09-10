% plotLRrunCorr.m
%
% Function for looking at correlation in activity between left and right
%  cells in calcium imaging.
% Prompts user to select pData file through GUI and outputs graph showing
%  two dF/F traces, FicTrac yaw and forward velocities, and running
%  correlation. Also shows overall correlation b/w two cells, during whole
%  trial and only during moving bouts
%
% INPUTS:
%   corrWin - length of window, in seconds, over which to compute running
%       correlation
%   prompts for pData file through GUI
%
% OUTPUTS:
%   none, but produces graph
%
% CREATED: 9/10/21 - HHY
%
% UPDATED:
%   9/10/21 - HHY
%
function plotLRrunCorr(corrWin)

    % prompt user to select pData file
    disp('Select pData file');
    [pDatFile, pDatPath] = uigetfile('*.mat', 'Select a pData file', ...
        pDataPath());
    
end
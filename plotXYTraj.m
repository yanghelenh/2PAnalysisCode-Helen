% plotXYTraj.m
%
% Function to plot fly's trajectory in x-y space, where time is color coded
%
% INPUTS:
%   xPos - x position
%   yPos - y position
%   t - time
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 12/12/18 HHY
% UPDATED: 12/12/18 HHY
%

function plotXYTraj(xPos, yPos,t)
    
    figure;
    
    scatter(xPos, yPos, 1, t);
    colorbar
    axis equal
end

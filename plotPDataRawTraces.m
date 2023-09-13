% plotPDataRawTraces.m
%
% Quick script to plot raw traces given user-selected pData file
% Written to plot traces for CSHL Neurobiology of Drosophila 2019 meeting
%  poster
%
% CREATED: 9/28/19 - HHY

%% plotting parameters
dFFscale = [-1 1];
fwdVelScale = [-10 40];
yawVelScale = [-600 600];
xScale = [100 130];


%% load data

% curDir = pwd;
% cd(pDataPath()); % go to pData folder
% 
% % ask user to select pData file
disp('Select a pData file to display.');
uPData = uigetfile('*pData.mat');
load(uPData);

%% generate plot
numSubplots = 4;


imgStartInd = find(img.t < xScale(1), 1, 'last');
imgEndInd = find(img.t > xScale(2), 1, 'first');

ftStartInd = find(fictrac.t < xScale(1), 1, 'last');
ftEndInd = find(fictrac.t > xScale(2), 1, 'first');

f = figure('Position', [0 0 1200 900]);

h{1} = subplot(numSubplots, 1, 1);
plot(img.t(imgStartInd:imgEndInd), img.filtDFF.left(imgStartInd:imgEndInd));
xlim(xScale);
ylim(dFFscale);
ylabel('Left dF/F');

h{2} = subplot(numSubplots, 1, 2);
plot(img.t(imgStartInd:imgEndInd), img.filtDFF.right(imgStartInd:imgEndInd));
xlim(xScale);
ylim(dFFscale);
ylabel('Right dF/F');

h{3} = subplot(numSubplots, 1, 3);
plot(fictrac.t(ftStartInd:ftEndInd), fictrac.fwdVel(ftStartInd:ftEndInd));
xlim(xScale);
ylim(fwdVelScale);
ylabel('Forward velocity (mm/s)');

h{4} = subplot(numSubplots, 1, 4);
plot(fictrac.t(ftStartInd:ftEndInd), fictrac.yawAngVel(ftStartInd:ftEndInd));
xlim(xScale);
ylim(yawVelScale);
ylabel('Yaw velocity (deg/s)');
xlabel('Time (s)');


% cd(curDir);
% plotCaVelExTraces.m
%
% Quick script to plot example traces calcium imaging traces (left and
%  right dF/F, yaw velocity, forward velocity). For making figure of
%  example traces for paper
%
% Adaptation of plotPDataRawTraces.m
%
% CREATED: 5/13/24 - HHY

%% plotting parameters
dFFscale = [-1 1];
fwdVelScale = [-20 40];
yawVelScale = [-750 750];
xScale = [175 195];


%% load data from pData file

% select pData file
[pDataFName, pDataDirPath] = uigetfile('*.mat', ...
    'Select pData file', pDataPath, 'MultiSelect', 'off');

pDataFullPath = [pDataDirPath filesep pDataFName];

load(pDataFullPath,'img','fictrac');

%% smooth FicTrac
% to match filtering of other data in paper
padLen = 100;
kernLen = 50;

yawSmo = gaussSmooth(fictrac.yawAngVel,padLen,kernLen);
fwdSmo = gaussSmooth(fictrac.fwdVel,padLen,kernLen);

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
plot(fictrac.t(ftStartInd:ftEndInd), yawSmo(ftStartInd:ftEndInd));
xlim(xScale);
ylim(yawVelScale);
ylabel('Yaw velocity (deg/s)');
xlabel('Time (s)');

h{4} = subplot(numSubplots, 1, 4);
plot(fictrac.t(ftStartInd:ftEndInd), fwdSmo(ftStartInd:ftEndInd));
xlim(xScale);
ylim(fwdVelScale);
ylabel('Forward velocity (mm/s)');

sgtitle(sprintf('%s',pDataFName));

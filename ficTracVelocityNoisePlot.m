% ficTracVelocityNoisePlot.m
%
% quick and dirty script to plot raw FicTrac yaw signal, showing where
% oscillations in velocity calculation come from (slight jitter in position
% signals are amplified)
%
% see FicTracVelocityNoise.fig figure saved in 2PAnalysis-Helen folder

cd('/Volumes/Samsung_2tb/Wilson Lab/RAW DATA/190113/fly01/fov01/trial01');

load('fictracDat.mat');

sampRate = 1/(median(diff(t)));

% moving average filter
avgWindow = 0.5; % in seconds
windowSize = avgWindow * sampRate;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smoYawAngVel = filtfilt(b,a, yawAngVel);

figure;
subplot(2,1,1);
plot(t, yawAngPos);
xlim([0,15]);
title('Unwrapped Yaw Position');
xlabel('Time (sec)');
ylabel('Yaw Position (deg)');

subplot(2,1,2);
plot(t,yawAngVel);
hold on;

plot(t,smoYawAngVel);

line([0,15],[0,0], 'Color', 'black'); % x-axis line

legend({'Angular Velocity', 'Moving Average Smoothed Angular Velocity'});

xlim([0,15]);
title('Yaw Angular Velocity');
xlabel('Time (sec)');
ylabel('Yaw Angular Velocity (deg/sec)');



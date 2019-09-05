% loadKernelYScaleYLabels.m
%
% Quick script to load y scales and y labels for kernel plotting, for when
% behavioral velocity units are all deg/s or for both mm/s and deg/s
%
% UPDATED: 9/4/19 - HHY

yScale = {[-4e-4 16e-4]...
    [-1e-5 2e-5]...
    [-1e-3 2e-3]...
    [-1e-5 4e-5]...
    [-0.5e-5 2.5e-5]...
    [-0.15 0.4]...
    [-20 20]...
    [-0.3 0.3]...
    [-2 10]...
    [-5 20]};

yLabels = {'dF/F unit per mm/s',...
    'dF/F unit per deg/s', ...
    'dF/F unit per mm/s', ...
    'dF/F unit per deg/s', ...
    'dF/F unit per deg/s', ...
    'mm/s per unit dF/F', ...
    'deg/s per unit dF/F', ...
    'mm/s per unit dF/F', ...
    'deg/s per unit dF/F', ...
    'deg/s per unit dF/F'};

degPerMM = 17.738631428198860;

yScaleAllDeg = {[-3e-5 9e-5]...
    [-1e-5 2e-5]...
    [-6e-5 10e-5]...
    [-1e-5 4e-5]...
    [-0.5e-5 2.5e-5]...
    [-3 10]...
    [-20 20]...
    [-5 5]...
    [-2 10]...
    [-5 20]};

yLabelsAllDeg = {'dF/F unit per deg/s',...
    'dF/F unit per deg/s', ...
    'dF/F unit per deg/s', ...
    'dF/F unit per deg/s', ...
    'dF/F unit per deg/s', ...
    'deg/s per unit dF/F', ...
    'deg/s per unit dF/F', ...
    'deg/s per unit dF/F', ...
    'deg/s per unit dF/F', ...
    'deg/s per unit dF/F'};
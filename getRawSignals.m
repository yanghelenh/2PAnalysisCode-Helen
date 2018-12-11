% getRawSignals.m
% 
% Function that computes average fluorescence over time in all the ROIs.
%  Outputs both average fluorescence with and without background
%  subtraction.
%  The functions drawROIs and selectBackground must have been run first to
%  define the ROIs and background region
%
% INPUT:
%   alSeries - aligned time series
%   roiMasks - masks for all ROIs
%
% OUTPUT:
%   avSignal - average pixel intensity in each ROI, vector length k masks
%   bksSignal - average ROI intensity minus average background intensity
%
% UPDATED: 12/10/18 HHY
%

function [avSignal, bksSignal] = getRawSignals(alSeries, roiMasks, bkMask)

    nFrames = size(alSeries,3); % number of frames in time series
    nMasks = length(roiMasks); % number of masks (i.e. num ROIs)
    
    % preallocate arrays for ROI time series
    avSignal = zeros(nMasks, nFrames);
    bksSignal = zeros(nMasks, nFrames);
    sizeROIMasks = zeros(nMasks,1);
    
    for k = 1:nMasks
        % number of pixels in each ROI
        sizeROIMasks(k) = sum(sum(roiMasks{k}));
    end

    % as above, but for background, not ROIs
    sizeBkMask = sum(sum(bkMask));

    % compute signal in each ROI
    for i = 1:nFrames
        A = double((alSeries(:,:,i))); % individual frame image
        bkMaskImg = A(bkMask); % background region

        for k = 1:nMasks
            maskImg = A(roiMasks{k}); % ROI region
            % average intensity in ROI (sum over space, divide by ROI size)
            avSignal(k,i) = sum(maskImg)./sizeROIMasks(k);
            % background subtraction (by average background intensity)
            bksSignal(k,i) = avSignal(k,i) - sum((bkMaskImg))./sizeBkMask;
        end
    end

end 
% computeWienerKernel.m
%
% Function to compute first-order Wiener kernel, as estimate of linear
%  filter, using the method described in French 1976, a la Clark et al.
%  2011.
% Deals with NaNs in input or output by excluding any segment that contains
%  any.
% Note: input and output must be same length and have same sampling rate.
%
% INPUT:
%   input - stimulus signal (i.e. signal whose autocorrelation appears in
%       the denominator of the kernel computation)
%   output - response signal (i.e. signal that is being cross-correlated in
%       the numerator of the kernel computation).
%   sampRate - sampling rate of input and output signals
%   winLen - length of window, in samples, that averages are computed over
%   lowPassCutoff - cutoff of lowpass filter to apply to numerator, use 0 
%       if no filtering desired
%
% OUTPUT:
%   kernel - time-domain estimate of kernel
%
% CREATED: 4/1/19
% UPDATED: 4/1/19
%

function kernel = computeWienerKernel(input, output, sampRate, winLen,...
    lowPassCutoff)

    % get overlapping segments of input and output of length window
    % total number of segments
    numSeg = length(input) - winLen + 1;
    % preallocate
    inputSeg = zeros(numSeg, winLen);
    outputSeg = inputSeg;
    
    % extract segments
    for i = 1:numSeg
        endInd = i + winLen - 1;
        inputSeg(i,:) = input(i:endInd);
        outputSeg(i,:) = output(i:endInd);
    end
    
    % remove any segments with NaN in either input or output
    % gets indicies of segments with at least 1 NaN
    inputNaNs = find(sum(isnan(inputSeg), 2)); 
    outputNaNs = find(sum(isnan(outputSeg), 2));
    % segments that have NaNs in input or output or both
    nanInd = union(inputNaNs, outputNaNs);
    
    % remove segments with NaNs
    inputSeg(nanInd,:) = [];
    outputSeg(nanInd,:) = [];
    
    % compute Wiener Kernel
    numSeg = size(inputSeg, 1); % update numSeg after NaN segments removed
    
    % get fft of each segment, input and output
%     n = (window-1) * 2; % new input length for zero-padding
    n = 2^nextpow2(winLen);
    inFFT = fft(inputSeg, n, 2);
    outFFT = fft(outputSeg, n, 2);
    inConj = conj(inFFT); % complex conjugate of input
    
    % first order Wiener Kernel is <Y(w)X*(w)>/<X(w)X*(w)>
    numerator = mean(outFFT .* inConj, 1);
    denominator = mean(inFFT .* inConj, 1);
    
    % low pass filter numerator
    if (lowPassCutoff > 0) % apply low pass filtering to numerator
        % actual frequencies for first half
        halfFreq = sampRate * (0:(n/2))/n; 
        % actual frequencies for whole fft
        freq = [halfFreq(1:(end-1)) fliplr(halfFreq(1:(end-1)))];
        % weights at each frequency, for first order filter with specified
        %  low pass cutoff
        filtWeights = 1 ./ (1 + (freq./(2 * pi * lowPassCutoff)));
        
        % low pass filter the numerator
        lpfNumerator = numerator .* filtWeights;
        
        % compute 1st order wiener kernel
        kernelFFT = lpfNumerator ./ denominator;
          
    else % no filtering
        % compute 1st order wiener kernel 
        kernelFFT = numerator ./ denominator;
    end
    
    % transform frequency domain wiener kernel into time domain
    kernel = ifft(kernelFFT,n);
    kernel = kernel(1:winLen); % remove zero padding
    
end
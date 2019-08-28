% plotFreqSpec.m
%
% Function that takes in time series and its sampling rate and plots the
%  one-sided frequency spectrum
%
% INPUTS:
%   in - time series to get frequency spectrum of
%   sampRate - sampling rate of time series
%
% OUTPUTS:
%   none, but generates frequency spectrum plot
%
% CREATED: 8/27/19 - HHY
% UPDATED: 8/27/19 - HHY

function plotFreqSpec(in, sampRate)

    % length of input
    lenIn = length(in);
    
    n = 2^nextpow2(lenIn);
    
    % fft of input
    fftIn = fft(in, n);
    
    % x-axis frequencies
    f = sampRate*(0:(n/2))/n;
    
    % only keep one half of frequency spectrum, magnitudes only
    P2 = abs(fftIn/n);
    P1 = P2(1:n/2+1);
    
    % plot
    figure;
    plot(f, P1);
    xlabel('Frequency (Hz)');

end
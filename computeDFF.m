% computeDFF.m
%
% Function that takes in background-subtracted fluorescence signal over
%  time, performs bleaching correction (sum of two exponentials), and
%  returns dF/F. 
% For bleaching correction, fits sum of two exponentials to portion of
%  fluorescence curve taken as input and uses it as F0 in dF/F = (F-F0)/F0
%  calculation. Which portion of fluorescence curve is used for fitting
%  effectively defines where dF/F=0 sits.
%
% INPUT
%   bksSignal - background-subtracted fluorescence signal over time, for
%       single ROI
%   imgFrameTimes - times at which bksSignal values were captured
%   baselineSignal - background-subtracted fluorescence signal over time,
%       for same ROI, but can be subset of bksSignal. Used for fitting
%       during bleaching correction.
%   baselineTimes - times at which baselineSignal values were captured
%
% OUTPUT
%   dFF - bleaching corrected dF/F trace for single ROI; as vector
%
function dFF = computeDFF(bksSignal, imgFrameTimes, baselineSignal, ...
    baselineTimes)    

    % make column vectors if needed
    if isrow(bksSignal)
        bksSignal = bksSignal';
    end 
    if isrow(baselineSignal)
        baselineSignal = baselineSignal';
    end 
    if isrow(imgFrameTimes)
        imgFrameTimes = imgFrameTimes';
    end
    if isrow(baselineTimes)
        baselineTimes = baselineTimes';
    end

    try
        fitobj = fit(baselineTimes, baselineSignal, fittype('exp2'));
    catch err 
        if(strcmp(err.identifier,'curvefit:fit:nanComputed') ||...
            strcmp(err.identifier,'curvefit:fit:xDataMustBeColumnVector'))
        else
            rethrow(err);
        end
    end 
    
    baselineF = fitobj(imgFrameTimes);
    dFF = bksSignal./baselineF - 1;  % equivalent to (F - F0)/F0
    
    clear fitobj

end 
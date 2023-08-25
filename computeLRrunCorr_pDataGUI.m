% computeLRrunCorr_pDataGUI.m
%
% Function to compute the running correlation between the left and right
%  cells' dF/F, given a window size to compute the correlation over
% Runs on all pData files selected through GUI
% Saves correlation and window size back into selected pData file
% Correlation computed on img.filtDFF values
%
% INPUTS:
%   corrWin - length of window, in seconds, over which to compute running
%       correlation
%
% OUTPUTS:
%   none, but saves back into original pData file
%
% CREATED:
%   12/13/22 - HHY
%
% UPDATED:
%   12/13/22 - HHY
%
function computeLRrunCorr_pDataGUI(corrWin)

    % prompt user to select pData files
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];
    
        % load data from pData file
        load(pDataFullPath, 'img', 'fictrac', 'name');
    
    
        % compute running correlation b/w left and right cells
        
        % get size of window, as odd number of frames
        ifi = median(diff(img.t)); % interframe interval
        corrWinFr = floor(corrWin / ifi);
        % if even, add 1
        if ~(mod(corrWinFr,2))
            corrWinFr = corrWinFr + 1;
        end
        halfWinSize = floor(corrWinFr / 2);
        
        
        % preallocate
        runCorr = zeros(size(img.filtDFF.left));
        
        % loop through all frames; window centered on that imaging frame
        for j = 1:length(img.filtDFF.left)
            startInd = j - halfWinSize; % get start index of window
            % window start would be before start of trial, just set to first
            %  frame
            if (startInd < 1)
                startInd = 1;
            end
            
            endInd = j + halfWinSize; % get end index of window
            % window start would be after end of trial, set to end
            if (endInd > length(img.filtDFF.left))
                endInd = length(img.filtDFF.left);
            end
            
            % dF/F snippet for left cell
            leftSeg = img.filtDFF.left(startInd:endInd);
            % dF/F snippet for right cell
            rightSeg = img.filtDFF.right(startInd:endInd);
            
            % get correlation b/w these two snippets
            thisCorrCoef = corrcoef(leftSeg,rightSeg);
            runCorr(j) = thisCorrCoef(1,2);   
        end
        
        % overall correlation
        overallCorr = corrcoef(img.filtDFF.left,img.filtDFF.right);
        overallCorr = overallCorr(1,2);

        % add runCorr into img.filtDFF
        img.filtDFF.corr = runCorr;

        % create struct for overallCorr and corrWin, in img struct
        img.corrParams.overallCorr = overallCorr;
        img.corrParams.corrWin = corrWin;

        % save updated img struct back into pData file
        save(pDataFullPath, 'img' ,'-append');
    end
end
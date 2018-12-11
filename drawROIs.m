% drawROIs.m
%
% Allows user to select ROIs from an image. Uses built-in Matlab function
% that lets you draw polygons on figures.
% 
% Inputs: 
%   avgIm - 2D image array (i.e. average image of time series)
%
% Outputs: 
%   roiMasks - cell array of logical array indicating pixels corresponding 
%     to ROIs
%

function roiMasks = drawROIs(avgIm)
 
    % display average image, full screen
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(avgIm,[],'InitialMagnification','fit'); 
    colormap('gray');
    
    title('Press ENTER when done selecting ROIs');
    
    % user selects ROIs
    done = 0;
    index = 1;
    while (~done)
        roiMasks{index} = roipoly; 
        done = waitforbuttonpress;
        index = index + 1;
    end 
    
    close(f);
   
end 


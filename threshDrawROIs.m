% threshDrawROIs.m
%
% Allows user to select ROIs from an image. First, user sets threshold in
%  dF/F; then, user draws rough region around ROI. ROI is intersection of
%  greater than threshold and within drawn region. Repeat for all ROIs.
%  Number of ROIs determined by user input at beginning.
% Uses built-in Matlab function that lets you draw polygons on figures.
% 
% Inputs: 
%   avgIm - 2D image array (i.e. average image of time series)
%
% Outputs: 
%   roiMasks - cell array of logical array indicating pixels corresponding 
%     to ROIs
%
% Created: 2/13/19 HHY
% Updated: 2/13/19 HHY


function roiMasks = threshDrawROIs(avgIm)

    % dF/F image, across space in average image
    % What does this mean? 
    %  (avg image)/(avg pixel intensity of avg image)-1 represents 
    %  deviation from the mean intensity. 
    devImg = avgIm./mean(avgIm(:))-1; 
    
    % threshold image at inital thresh
    thresh = 0;
    md = imdilate(devImg > thresh,ones(3));
    mdr = imerode(md,ones(3));
    
    % display figure
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % average image
    subplot(2, 2, 1);
    imshow(avgIm, [], 'InitialMagnification','fit'); 
    colormap('gray');
    title('Average image');
    
    % deviation from mean intensity image
    subplot(2, 2, 2);
    imshow(devImg, [], 'InitialMagnification','fit'); 
    colormap('gray');
    title('Deviation from mean image');
    
    % histogram of intensity deviations
    subplot(2, 2, 4);
    histogram(devImg(:),200); 
    set(gca,'yscale','log');
    title('Histogram of deviations from mean intensity');
    
    % average image (will plot thresholded mask on top of)
    subplot(2, 2, 3);
    imshow(avgIm, [], 'InitialMagnification','fit'); 
    colormap('gray');
    title('Average image');  
    
    % click on figure; goes back to command line and request for num ROIs
    waitforbuttonpress;
    commandwindow
    
    % ask for number of ROIs from user
    numROIs = [];
    prompt = 'Number of ROIs: ';

    % repeat until number is entered
    while (isempty(numROIs))
        ui = input(prompt,'s');
        numROIs = str2num(ui);
    end
        
    % color mask 
    xPixels = size(avgIm, 1);
    yPixels = size(avgIm, 2);
    rCMask = zeros(xPixels, yPixels, 3);
    cm = colormap(lines);
    
    % loop through all ROIs
    for i = 1:numROIs

        % initialize loop parameters
        DONE = 0;
        thresh = 0;
        but = 1;
            
        % establish mask color for this ROI
        if (i>64)
            cmInd = mod(i,65)+1;
        else
            cmInd = i;
        end
        curColor = cm(cmInd,:);
        
        % initialize mask
        CMask = zeros(xPixels, yPixels, 3);
            
        while ~DONE
            
            threshkeep = thresh; % just keep it, in case...
            buttonkeep = but;
            md = imdilate(devImg > thresh,ones(3));
            mdr = imerode(md,ones(3));
            
            subplot(2, 2, 3);
            imshow(avgIm, [], 'InitialMagnification','fit'); 
            colormap('gray');
            title('Average image');

            % get color mask for thresholded cells
            CMask = cat(3, curColor(1).*mdr, ...
                curColor(2).*mdr, curColor(3).*mdr);
            
            % plot color mask
            hold on;
            h = imshow(CMask);
            set(h,'AlphaData',0.5);
            
            subplot(2,2,4);
            [thresh,~,but] = ginput(1); % select
            DONE = isempty(thresh); % exit loop when return pressed...
        end
        
        % mask image just by thresholding
        threshMask = devImg > threshkeep;
        
        % allow user to draw mask
        subplot(2, 2, 3);
        drawnMask = roipoly;
        
        roiMasks{i} = threshMask .* drawnMask;
        roiMasks{i} = logical(roiMasks{i});
        
        % plot ROI mask
        curMask = cat(3, curColor(1).*roiMasks{i}, ...
            curColor(2).*roiMasks{i}, curColor(3).*roiMasks{i});
        rCMask = rCMask + curMask;
        
        subplot(2, 2, 3);
        imshow(avgIm, [], 'InitialMagnification','fit'); 
        colormap('gray');
        title('Average image');
        hold on;
        h = imshow(rCMask);
        set(h,'AlphaData',0.5);
        
        waitforbuttonpress;
        
    end
    
    close(f);
   
end 


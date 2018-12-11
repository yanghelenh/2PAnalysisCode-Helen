% selectBackground.m
%
% Function that allows the user to select the background region from the 
%  time series's average image. Displays a histogram of pixel intensities
%  (normalized as deviation from the mean) and allows the user to select a
%  static threshold. Background region is then used for background
%  subtraction of fluorescence values before dF/F is computed.
%
% INPUT:
%   AV - the time series's average image, as 2D array
%
% OUTPUT:
%   bkMask - logical array indicating background regions
%

function [bkMask] = selectBackground(AV)

    % select background region    

    % dF/F image, across space in average image
    % What does this mean? 
    %  (avg image)/(avg pixel intensity of avg image)-1 represents 
    %  deviation from the mean intensity. 
    AV = double(AV);
    D = AV./mean(AV(:))-1; 

    DONE = 0;
    thresh = 0;
    but = 1;
    md = imdilate(D > thresh,ones(3));
    mdr = imerode(md,ones(3));
    
    f = figure; 
    f.Name = 'Set threshold for Background. Press enter when done.';
    
    while ~DONE
        threshkeep = thresh; % just keep it, in case...
        buttonkeep = but;
        md = imdilate(D > thresh,ones(3));
        mdr = imerode(md,ones(3));
        
        % average image
        subplot(3,1,3); 
        imagesc(AV); colormap gray;
        title('average image');
        
        % histograms of deviations from mean intensity
        subplot(3,1,1);
%         hist(D(:), 200);
        histogram(D(:),200); 
        set(gca,'yscale','log');
        title('histogram of deviations from mean intensity');
        
        % background vs. cell boundary
        subplot(3,1,2);
        imagesc(mdr); %colormap('gray');
        title('background vs. cell boundary');
        subplot(3,1,1);
        [thresh,dum,but] = ginput(1); % select
        DONE = isempty(thresh); % exit loop when return pressed...
        
    end
    close(f);
    
    bkMask = D < threshkeep; % background mask binary image

end 
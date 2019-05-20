% createImgLegFictracMovie.m
%
% script to generate movie showing dF/F, legVid, FicTrac variables

trialPath = ...
    '/Volumes/Samsung_2tb/Wilson Lab/RAW DATA/190507/fly01/fov02/trial01';

tVid = [520 540]; % start and end times, in seconds

cd(trialPath);

% load relevant data
load('imDat.mat', 'alignedSeries', 'frameStartTimes');
load('pData.mat', 'dFFs');
load('legVidDat.mat', 'legVidFrameTimes');
load('fictracDat.mat', 't', 'fwdVel', 'yawAngVel');
ficTracTimes = t;

% process imaging data
% preallocate
% imgDFFs = zeros(size(alignedSeries.ch1));
% for i = 1:size(imgDFFs, 1)
%     for j = 1:size(imgDFFs, 2)
%         ch1Trace = double(squeeze(alignedSeries.ch1(i,j,:)));
%         ch1dFFTemp = computeDFF(ch1Trace, ...
%             frameStartTimes, ch1Trace, ...
%             frameStartTimes) + 1;
%         ch2Trace = double(squeeze(alignedSeries.ch2(i,j,:)));
%         ch2dFFTemp = computeDFF(ch2Trace, ...
%             frameStartTimes, ch2Trace, ...
%             frameStartTimes) + 1;
%         ch1DivCh2Temp = ch1dFFTemp ./ ch2dFFTemp - 1; 
%         
%         imgDFFs(i,j,:) = ch1DivCh2Temp;    
%     end
%     fprintf('dF/F calculation on row %d \n', i);
% end

% smooth alignedSeries first
smAlignedSeries.ch1 = zeros(size(alignedSeries.ch1));
smAlignedSeries.ch2 = zeros(size(alignedSeries.ch2));
smoWindow = 4;

for l = 1:size(alignedSeries.ch1,3)
    if (l <= smoWindow/2)
        smAlignedSeries.ch1(:,:,l) = alignedSeries.ch1(:,:,l);
        smAlignedSeries.ch2(:,:,l) = alignedSeries.ch2(:,:,l);
    elseif (l >= size(alignedSeries.ch1,3)-smoWindow/2)
        smAlignedSeries.ch1(:,:,l) = alignedSeries.ch1(:,:,l);
        smAlignedSeries.ch2(:,:,l) = alignedSeries.ch2(:,:,l);
    else
        smAlignedSeries.ch1(:,:,l) = ...
            mean(alignedSeries.ch1(:,:,(l-smoWindow/2):(l+smoWindow/2)),3);
        smAlignedSeries.ch2(:,:,l) = ...
            mean(alignedSeries.ch2(:,:,(l-smoWindow/2):(l+smoWindow/2)),3);
    end
end

% calculate dF/F on smoothed aligned series
% smImgDFFs = zeros(size(alignedSeries.ch1));
% for i = 1:size(smImgDFFs, 1)
%     fprintf('dF/F calculation on row %d \n', i);
%     for j = 1:size(smImgDFFs, 2)
%         ch1Trace = squeeze(smAlignedSeries.ch1(i,j,:));
%         ch1dFFTemp = computeDFF(ch1Trace, ...
%             frameStartTimes, ch1Trace, ...
%             frameStartTimes) + 1;
%         ch2Trace = squeeze(smAlignedSeries.ch2(i,j,:));
%         ch2dFFTemp = computeDFF(ch2Trace, ...
%             frameStartTimes, ch2Trace, ...
%             frameStartTimes) + 1;
%         ch1DivCh2Temp = ch1dFFTemp ./ ch2dFFTemp - 1; 
%         
%         smImgDFFs(i,j,:) = ch1DivCh2Temp;    
%     end
% end

imgStartInd = find(frameStartTimes >= tVid(1), 1, 'first');
imgEndInd = find(frameStartTimes < tVid(2), 1, 'last');

% smImgDFFs = zeros(size(imgDFFs,1), size(imgDFFs,2), ...
%     imgEndInd - imgStartInd + 1);
smoWindow = 6; 

% for l = imgStartInd:1:imgEndInd
%     smImgDFFs(:,:,l-imgStartInd+1) = ...
%         mean(imgDFFs(:,:,(l-smoWindow/2):(l+smoWindow/2)),3);
% end

smRaw.ch1 = zeros(size(alignedSeries.ch1,1), size(alignedSeries.ch1,2), ...
    imgEndInd - imgStartInd + 1);
smRaw.ch2 = zeros(size(alignedSeries.ch2,1), size(alignedSeries.ch2,2), ...
    imgEndInd - imgStartInd + 1);

for l = imgStartInd:1:imgEndInd
    smRaw.ch1(:,:,l-imgStartInd+1) = ...
        mean(alignedSeries.ch1(:,:,(l-smoWindow/2):(l+smoWindow/2)),3);
    smRaw.ch2(:,:,l-imgStartInd+1) = ...
        mean(alignedSeries.ch2(:,:,(l-smoWindow/2):(l+smoWindow/2)),3);  
    smRaw.ch1(:,:,l-imgStartInd+1) = ...
        imgaussfilt(smRaw.ch1(:,:,l-imgStartInd+1),1);
    smRaw.ch2(:,:,l-imgStartInd+1) = ...
        imgaussfilt(smRaw.ch2(:,:,l-imgStartInd+1),1);
end
% smRaw.norm = (smRaw.ch1 ./ mean(smRaw.ch1,3)) ./ ...
%     (smRaw.ch2 ./ mean(smRaw.ch2,3));
smRaw.norm = smRaw.ch1 ./ smRaw.ch2;

% plot smooth raw norm
% minInt = squeeze(min(min(min(smRaw.norm))));
% maxInt = squeeze(max(max(max(smRaw.norm))));
maxImgInt = 2800;
minImgInt = 280;

% imaging frame times
imgFrameTimes = frameStartTimes(imgStartInd:imgEndInd) - ...
    frameStartTimes(imgStartInd);

% img dF/F scale
dFFScale = [-0.5 0.5];


% FicTrac variables
sampRate = 1/(median(diff(ficTracTimes)));
avgWindow = 0.3;
smoYawAngVel = moveAvgFilt(yawAngVel, sampRate, avgWindow);
smoFwdVel = moveAvgFilt(fwdVel, sampRate, avgWindow);


% legVid variables
legStartInd = find(legVidFrameTimes >= tVid(1), 1, 'first');
legEndInd = find(legVidFrameTimes < tVid(2), 1, 'last');
legVidPath = [pwd filesep 'rawLegVid'];
legVidImgSize = [416 416];

legVidImg = zeros(legVidImgSize(1),legVidImgSize(2), ...
    legEndInd - legStartInd + 1);
for i = legStartInd:1:legEndInd
    legVidFile = sprintf('legVid-%i.tiff', i-1);
    legVidImg(:,:,i - legStartInd + 1) = imread(...
        [legVidPath filesep legVidFile]);
end
legVidMinInt = 50;
legVidMaxInt = 150;

legVidFrameRate = 1/(median(diff(legVidFrameTimes)));
imgFrameRate = 1/(median(diff(frameStartTimes)));

% get indicies to display on each frame
imgDispInd = zeros(1,legEndInd - legStartInd + 1);
ficTracDispInd = zeros(1,legEndInd - legStartInd + 1);
for i = 1:(legEndInd - legStartInd + 1)
    imgDispInd(i) = find(frameStartTimes >= ...
        legVidFrameTimes(i + legStartInd - 1), 1, 'first') - ...
        imgStartInd + 1;
    ficTracDispInd(i) = find(ficTracTimes == ...
            legVidFrameTimes(i + legStartInd - 1), 1, 'first');  
end

ficTracStartInd = ficTracDispInd(1);
ficTracFigTimes = ficTracTimes(ficTracDispInd) - ficTracTimes(ficTracStartInd);
smoFwdVelFig = smoFwdVel(ficTracDispInd);
smoAngVelFig = smoYawAngVel(ficTracDispInd);

fwdVelScale = [-10 20];
angVelScale = [-500 500];

timeScale = [0 20];

% make figure
figFrames(legEndInd - legStartInd + 1) = struct('cdata',[],'colormap',[]);
% figFrames(200) = struct('cdata',[],'colormap',[]);
f = figure;
set(f,'Position',[10,10,1000,700]);
for i = 1:(legEndInd - legStartInd + 1)
    imgInd = imgDispInd(i);
    ficTracInd = ficTracDispInd(i);
    
    % calcium imaging
    imgAx = subplot(4, 3, [1 4]);
    imagesc(flipud(smRaw.ch1(:,:,imgInd)),[minImgInt,maxImgInt]); 
    axis equal;
    axis tight;
    colormap(imgAx,'parula');
    set(gca,'XTick',[],'YTick',[]);
    imgPos = get(imgAx,'Position');
    newImgPos = [imgPos(1)-0.06, imgPos(2)-0.06, imgPos(3)*1.5, imgPos(4)*1.5];
    set(imgAx,'Position',newImgPos);
    
    
    % leg video
    legAx = subplot(4, 3, [7, 10]);
    imagesc(legVidImg(:,:,i),[legVidMinInt legVidMaxInt]);
    axis equal;
    axis tight;
    colormap(legAx,'gray');
    set(gca,'XTick',[],'YTick',[]);
    legPos = get(legAx,'Position');
    newLegPos = [legPos(1)-0.06, legPos(2)-0.06, legPos(3)*1.5, legPos(4)*1.5];
    set(legAx,'Position',newLegPos);
        
    % calcium imaging traces, left
    leftDFFax = subplot(4,3,[2,3]);
    plot(imgFrameTimes(1:imgInd), ...
        dFFs(1,imgStartInd:(imgStartInd + imgInd - 1)), 'b', ...
        'LineWidth',1.5);
    ylim(dFFScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('dF/F');
    leftPos = get(leftDFFax,'Position');
    newLeftPos = [leftPos(1)+0.06, leftPos(2)+0.06, leftPos(3), leftPos(4)*0.9];
    set(leftDFFax,'Position',newLeftPos);
    title('Left Cell');
    
    % calcium imaging traces, right
    rightDFFax = subplot(4,3,[5,6]);
    plot(imgFrameTimes(1:imgInd), ...
        dFFs(2,imgStartInd:(imgStartInd + imgInd - 1)), 'r', ...
        'LineWidth',1.5);
    ylim(dFFScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('dF/F');
    rightPos = get(rightDFFax,'Position');
    newRightPos = [rightPos(1)+0.06, rightPos(2)+0.06, rightPos(3), rightPos(4)*0.9];
    set(rightDFFax,'Position',newRightPos);   
    title('Right Cell')
    
    % FicTrac forward velocity
    fwdVelAx = subplot(4,3,[8,9]);
    plot(ficTracFigTimes(1:i),...
        smoFwdVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(fwdVelScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('mm/s');
    fwdVelPos = get(fwdVelAx,'Position');
    newFVPos = [fwdVelPos(1)+0.06, fwdVelPos(2), fwdVelPos(3), fwdVelPos(4)*0.9];
    set(fwdVelAx,'Position',newFVPos);   
    title('Forward Velocity');
    
    % FicTrac yaw velocity
    yawAx = subplot(4,3,[11,12]);
    plot(ficTracFigTimes(1:i),...
        smoAngVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(angVelScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('deg/s');
    yawPos = get(yawAx,'Position');
    newYawPos = [yawPos(1)+0.06, yawPos(2), yawPos(3), yawPos(4)*0.9];
    set(yawAx,'Position',newYawPos);   
    title('Rotational Velocity');   
    
    drawnow;
    figFrames(i) = getframe(gcf); 
end

% save video
vidName = [pwd filesep 'figureVid.mp4'];

v = VideoWriter(vidName, 'MPEG-4');
v.FrameRate = legVidFrameRate;
open(v);
writeVideo(v, figFrames);
close(v);


% save('vidVars.mat', 'imgDFFs','-v7.3');




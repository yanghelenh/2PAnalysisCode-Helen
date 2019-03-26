% selectDroppedFicTrac.m
%
% Function to manually select times in experiment where FicTrac was not
%  tracking the ball. Save those indicies into fictracDat.mat, to allow use
%  in future analyses.
%
% INPUT: 
%   none, but user selects trial folder to process
%
% OUTPUT:
%   none, but saves dropInd variable into fictracDat.mat
%
% Created: 3/26/19 HHY
% Updated: 3/26/19 HHY
%

function selectDroppedFicTrac()

    T_RANGE = 50; % how many seconds of FicTrac data to display at once

    % ask user to select trial folder
    disp('Select a trial folder analyze.');
    uTrialPath = uigetdir;
    curDir = pwd;
    cd(uTrialPath)
    
    fprintf('Selecting FicTrac dropped times in %s \n', uTrialPath);
    
    % check that selected folder has fictracDat.mat file
    trialPathFiles = dir(uTrialPath);
    trialPathFileNames = extractfield(trialPathFiles, 'name');
    hasFictracDat = sum(strcmp(trialPathFileNames, 'fictracDat.mat'));
    
    if (hasFictracDat)
        % get variables in fictracDat
        ftDatStrct = whos('-file', 'fictracDat.mat');
        ftDatVars = extractfield(ftDatStrct, 'name');
        
        % check whether dropped frames previously selected; if so, prompt
        %  user to ask whether they want to overwrite
        ftDropSelected = sum(strcmp(ftDatVars, 'dropInd'));
        
        if (ftDropSelected) 
            prompt = ['FicTrac dropping times have already been '...
                'selected for this trial.' ' Overwrite? (Y/N) '];
            ui = input(prompt,'s');
            if (~strcmpi(ui, 'Y')) 
                % stop running this function. don't overwrite 
                disp('Ending selectDroppedFicTrac. Nothing overwritten');
                cd(curDir); % return to previous directory
                return;
            end
        end
        
        % load fictracDat.mat, only needed variables
        load('fictracDat.mat', 'fwdVel', 'yawAngVel', 't');
        
        % pre-generate dropInd variable and save into fictracDat
        dropInd = [];
        save([uTrialPath filesep 'fictracDat.mat'], ...
            'dropInd', '-append');
        
        % generate plot for selecting times of FicTrac dropping
        
        % determine starting xLim on all plots
        tLims = [0 T_RANGE];
        
        % end time of fictrac data
        xMax = t(end);

        f = figure;
        hold on;
        
        % plot fwd velocity
        fwdVelSubplot = subplot(2, 1, 1);
        fwdVelPlot = plot(t, fwdVel);
        xlim(tLims);
%         ylim([-10 30]);
        xlabel('Time (sec)');
        ylabel('mm/sec');
        title('Forward Velocity');
        
        % plot yaw angular velocity
        yawVelSubplot = subplot(2, 1, 2);
        yawVelPlot = plot(t, yawAngVel);
%         ylim([-500 500]);
        xlim(tLims);
        xlabel('Time (sec)');
        ylabel('deg/sec');
        title('Rotational Velocity');
        
        % link x-axes
        linkaxes([fwdVelSubplot, yawVelSubplot], 'x');
        
        % create brush object
        brushobj = brush(f);
        brushobj.Enable = 'off'; % start with brushobj off
        
        % slider for scrolling x-axis
        tSlider = uicontrol(f, 'Style', 'slider', 'Position', ...
            [20 10 400 20]);
        tSlider.Value = 0;
        tSlider.Callback = @updateTLim;

        % push button to activate/inactivate brushing to select dropped
        %  times
        brushBut = uicontrol(f, 'Style', 'togglebutton', 'Position', ...
            [500 10 50 30]);
        brushBut.Callback = @brushButtonPush;
        brushBut.Value = brushBut.Min; % start with button deselected
        
        cd(curDir)

        
        
    else
        disp('Selected trial folder does not contain fictracDat.mat file');
        cd(curDir);
        return;
    end

    % functions for interacting with plot
    % function to change time axis with slider
    function updateTLim(src, event)
        xlim(fwdVelSubplot, [tSlider.Value * (xMax-T_RANGE),...
                tSlider.Value * (xMax-T_RANGE) + T_RANGE]);
        xlim(yawVelSubplot, [tSlider.Value * (xMax-T_RANGE),...
                tSlider.Value * (xMax-T_RANGE) + T_RANGE]); 
    end

    % function to start/stop brushing, save values
    function brushButtonPush(src, event)
        % if toggle button is selected
        if (brushBut.Value == brushBut.Max)
            % enable brushing
            brushobj.Enable = 'on';
        
        % if toggle button is deselected
        elseif (brushBut.Value == brushBut.Min)
            % disable brushing
            brushobj.Enable = 'off';
            
            % get brushed values
            brushFwdVelInd = find(get(fwdVelPlot, 'BrushData'));
            brushYawVelInd = find(get(yawVelPlot, 'BrushData'));
            
            % FicTrac dropped values are union of values selected on 2
            %  plots
            dropInd = union(brushFwdVelInd, brushYawVelInd);
            
            % save into fictracDat.mat
            save([uTrialPath filesep 'fictracDat.mat'], ...
                'dropInd', '-append');
        end
    end
    
end
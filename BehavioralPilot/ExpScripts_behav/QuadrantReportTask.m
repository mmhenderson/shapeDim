% QUADRANT REPORT TASK - SHAPE SPACE 
% Subject sees a sequence of images - must report the categorical
% "quadrant" of each presented image. 
% Each quadrant has a corresponding prototype which is used to remind the
% subject of the relevant quadrants. 

% WRITTEN FOR BEHAVIOR ROOMS (B,C,D)

% parameters
    % Image size 16 degrees (fills almost entire screen)
    % Image onscreen for 1 second
    % Response period for 1 second (blank)
    % ITI 4 seconds (incl blank period)
    
    % 36 total images presented
    
    % total time = 2s blank start + 36*2*(1+1+2) = 290 seconds
    % 290/60 = 4.8 minutes
    
try
    
    echo off
    clear 
    close all hidden
       
    scannerlaptop = 0;
    saveremote = 1;
    
    expdir = pwd;
    filesepinds = find(expdir==filesep);
    root = expdir(1:filesepinds(end)-1);
    
    % set up paths for saving
    if scannerlaptop || IsWindows
        % if scanner laptop or debugging, save to a subfolder here.
        datadir_local = fullfile(root, 'Data');
        datadir_remote = datadir_local;
    else
        % linux, behavior rooms - save to a different folder.
        datadir_local = '/home/pclexp/Documents/Maggie/shapeDim/Data/';
        datadir_remote = '/mnt/pclexp/Maggie/shapeDim/Data/';
    end
    
    % where will i look for images?
    p.imagedir=fullfile(root, 'Stimuli/AmpGrid1/');
    
    %% Collect information about the subject, the date, etc.
    p.Subject = input('Initials? (default is tmp)--> ','s'); 
    if isempty(p.Subject)
        p.Subject = 'tmp'; 
    end 
    
    SubNum = input('Subject number?  --> ','s'); 
    if isempty(SubNum)
        error('enter a valid subject number')
    end
    p.SubNum = sprintf('%02d', str2double(SubNum));
    
    isTrn = input('Training run? (default is 0)  --> ','s'); 
    if isempty(isTrn)
        isTrn = '0';
    end
    p.Training = str2double(isTrn);
    
    RunNum = input('Run?  --> ','s');
    if isempty(RunNum)
        error('enter a valid run number')
    end
    p.runNumGlobal = str2double(RunNum);
    
    % difficulty: 1=easiest, 13=hardest
    Diff = input('Difficulty? (1-13, default is 6) --> ','s');
    if isempty(Diff)
        p.difficulty=6;
    else
        p.difficulty = str2double(Diff);
        if mod(p.difficulty,1) || p.difficulty<1 || p.difficulty>13
            error('difficulty must be an integer between 1 and 13')
        end
    end
    
    rng('default')
    t.MySeed = sum(100*clock); 
    rng(t.MySeed);
    
    if p.Training
        p.expName = 'ShapeQuadTask_TRAINING';
    else
        p.expName = 'ShapeQuadTask';
    end
    
    %% initialize my data file

    % save a copy of the currently running script here (in text format) so
    % we can look at it later if issues arise
    p.MyScript = fileread([mfilename('fullpath'),'.m']);
    
    % make sure my data dir exists
    if ~exist(datadir_local,'dir')
        mkdir(datadir_local)
    end
    if saveremote && ~exist(datadir_remote,'dir')
        mkdir(datadir_remote)
    end
    t.TheDate = datestr(now,'yymmdd'); %Collect todays date (in t.)
    t.TimeStamp = datestr(now,'HHMM'); %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
    
   p.fnsave_local = fullfile(datadir_local, ['S', p.SubNum, '_' p.expName '_' t.TheDate '.mat']);
    p.fnsave_remote = fullfile(datadir_remote, ['S', p.SubNum, '_' p.expName '_' t.TheDate '.mat']);
    
    if exist(p.fnsave_local,'file')
        load(p.fnsave_local);
        if p.runNumGlobal~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.runNumGlobal~=1            
        error('No data exists yet for this subject, check your run number')
    end
    
    p.rndseed = round(sum(100*clock));
    rng(p.rndseed);
    
    if p.Training
%         p.nTrials = 1;
        p.nTrials = 36;
    else
%         p.nTrials = 36*2;
        p.nTrials = 36;
%         p.nTrials = 1;
    end
    
    p.nIms = p.nTrials;
   
   %% set up the grid 
   
   % first I'm defining all possible images that we can use in this task.
    start = 0.2;    % min value along each axis
    stop = 4.8; % max value along each axis  
    step = 0.1;
    center = (stop-start)./2+start;
      
    all_pts = start:step:stop;  
    [gridx,gridy] = meshgrid(all_pts,all_pts);
    all_grid_points = [gridx(:),gridy(:)];
    
    % now taking out images at exactly the prototype locations so that we
    % never use these during task
    nsteps_main = 4;
    main_pts = round(linspace(start,stop, nsteps_main),1);
    proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
    proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
    proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
    assert(numel(proto_inds)==4)
    all_grid_points(proto_inds,:) = [];
  
    % Also taking out any images along quadrant boundaries because these
    % are ambiguous 
    bound_inds = find(any(all_grid_points==center,2));
    all_grid_points(bound_inds,:) = [];
    
    % Next, for each point in the full grid, define the difficulty level 
    % based on distance from the nearest quadrant boundary
    dist_from_bound = round(min(abs((all_grid_points-repmat(center, size(all_grid_points,1),2))),[],2),1);
     
    % bin these values into 13 "levels"
    [undist, ~ , dist_groups] = unique(round(dist_from_bound,1)); 
    
    % define the start of each bin
    bin_edges = round(fliplr([0.1:0.1:0.9,1.2:0.3:2.1]),1); 
    dist_groups_binned = zeros(size(dist_groups));
    for bb=1:numel(bin_edges)
        if bb>1
            inds = dist_from_bound>=bin_edges(bb) & dist_from_bound<bin_edges(bb-1);
        else
            inds = dist_from_bound>=bin_edges(bb);
        end
        dist_groups_binned(inds) = bb;
    end  
  
    
    % add a zero here to make the numbering work
%     undist = [0; undist];    
%     undist_binned = reshape(undist, 3, 8)';
%     dist_groups_binned = zeros(size(dist_groups));
%     for bb=1:size(undist_binned,1)
%         inds = ismember(round(dist_from_bound,1), undist_binned(bb,:));
%         dist_groups_binned(inds) = bb;
%     end    

    nperbin = sum(repmat(dist_groups_binned, 1, size(bin_edges,2))==repmat(1:size(bin_edges,2),size(dist_groups_binned,1),1),1);
    
    %% Now choose the images that we want to show on this block, randomly sampling from the correct bin.
    thisbin = find(dist_groups_binned==p.difficulty);    
    if numel(thisbin)<p.nIms
        % if this is the smallest (easiest) bin, we need to repeat images
        myinds = thisbin;
        myinds = [myinds; datasample(thisbin, p.nIms-numel(thisbin),'replace',false)];
    else
        % if it's a bigger bin, just take nIms images with no repeats
        myinds = datasample(thisbin, p.nIms, 'replace',false);
    end
    p.points = all_grid_points(myinds,:);

    % check this math
    dist_from_bound = round(min(abs((p.points-repmat(center, size(p.points,1),2))),[],2),1);
    if p.difficulty>1
        assert(all(dist_from_bound>=bin_edges(p.difficulty) & dist_from_bound<bin_edges(p.difficulty-1)))
    else
        assert(all(dist_from_bound>=bin_edges(p.difficulty)))
    end
        
    %% plot these if you want to see how the grid looks
%     close all
%     
%     figure;hold all;
%    
%     cols = parula(numel(bin_edges));
%     for dd=1:numel(bin_edges)
%         plot(all_grid_points(dist_groups_binned==dd,1), all_grid_points(dist_groups_binned==dd,2),...
%             'LineStyle','none','Marker','.',...
%             'MarkerFaceColor',cols(dd,:),'MarkerEdgeColor',cols(dd,:))        
%     end
%     axis equal
%     plot(p.points(:,1),p.points(:,2),'ro');
%     line([center,center],get(gca,'YLim'),'Color','k')
%     line(get(gca,'XLim'), [center,center], 'Color','k');
%     
%     set(gcf,'Color','w');
%     title(sprintf('Shape locations used for this block\n(Difficulty=%d)',p.difficulty));
    
    %% Make a randomized image sequence
    
    p.imlist = (1:p.nIms)';
    p.imlist = p.imlist(randperm(length(p.imlist)));
    p.imcoords = p.points(p.imlist,:);
    
    %define quadrant labels for each image (going counter-clockwise from top right)
    p.quadrant = zeros(size(p.imlist));
    p.quadrant(p.imcoords(:,1)>center & p.imcoords(:,2)>center) = 1;
    p.quadrant(p.imcoords(:,1)<center & p.imcoords(:,2)>center) = 2;
    p.quadrant(p.imcoords(:,1)<center & p.imcoords(:,2)<center) = 3;
    p.quadrant(p.imcoords(:,1)>center & p.imcoords(:,2)<center) = 4;

    %% set up my screen 
    
    InitializeMatlabOpenGL;  
    PsychImaging('PrepareConfiguration');
    Screen('Preference', 'SkipSyncTests', 0);
    AssertOpenGL; % bail if current version of PTB does not use
    PsychJavaTrouble;
    
    p.windowed = 0;
    s=max(Screen('Screens'));
    p.black = BlackIndex(s);
    p.white = WhiteIndex(s);
    % Open a screen
    Screen('Preference','VBLTimestampingMode',-1);  % for the moment, must disable high-precision timer on Win apps
    multiSample=0;
    if p.windowed
        Screen('Preference', 'SkipSyncTests', 1);
        [w, p.sRect]=Screen('OpenWindow', s, p.backColor,[50,50,800,600],[],[],multiSample);
    else

    % don't change this, matches the images!
    p.backColor = 0;

    [w, p.sRect]=Screen('OpenWindow', s, p.backColor,[],[],[],multiSample);
        HideCursor;
    end
    disp(p.sRect)
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % test the refresh properties of the display
    p.fps=Screen('FrameRate',w);          % frames per second
    p.ifi=Screen('GetFlipInterval', w);   % inter-frame-time
    if p.fps==0                           % if fps does not register, then set the fps based on ifi
        p.fps=1/p.ifi;
    end

    if scannerlaptop
        p.refreshRate = 60;
        % INNER BORE screen
        p.vDistCM = 47;
        p.screenHeightCM = 18;
%         % STAND UP SCREEN
%         p.screenHeightCM = 90; % in cm
%         p.vDistCM = 370; % in cm
    else
        
        if IsWindows
            % we're on windows lab PC, testing code
            p.refreshRate = 59;
        else
            % we're on linux machine in behav room
            p.refreshRate = 85;
        end
        % Behavior rooms (should all be identical)
        p.vDistCM = 49;
        p.screenHeightCM = 29;  
    end
    p.VisAngle = (2*atan2(p.screenHeightCM/2, p.vDistCM))*(180/pi); % visual angle of the whole screen
    
    % make sure the refreshrate is ok
    if abs(p.fps-p.refreshRate)>5
        Screen('CloseAll');
        disp('CHANGE YOUR REFRESH RATE')
        ListenChar(0);
        %clear all;
        return;
    end

    % if running the real experiment (not debugging), then hide cursor and set
    % priority high
    if ~p.windowed
        HideCursor; % Hide the mouse cursor
        % set the priority up way high to discourage interruptions
        Priority(MaxPriority(w));
    end

    %get center of screen
    p.centerPix = [(p.sRect(3) - p.sRect(1))/2, (p.sRect(4) - p.sRect(2))/2];
    p.fixSizeDeg = .2;
    p.fixColor = [0.8,0.8,0.8]*255;
%     p.fixColor = [0,0,0];
    p.stimHeightDeg = 14;
    p.protoStimHeightDeg = 4;
    p.protoDistDeg = 5;
    p.rectColor = [0.1,0.8,0.1]*255;
    p.rectWidthDeg = 0.1;
    % convert from degrees to pixel units
    p = deg2pix(p);
    p.fixSizePix = ceil(p.fixSizePix);
    %% Load the images

    for ii=1:p.nIms
        
        imfn = fullfile(p.imagedir, sprintf('Shape_%.2f_%.2f.png', p.points(ii,1),p.points(ii,2)));
        if exist(imfn,'file')
            im=imread(imfn);
        else
            error('image file %s not found!',imfn)
        end        
        
        allims(ii).name=imfn;
        allims(ii).imtext=Screen('MakeTexture',w,im);
        p.imfns{ii} = imfn;

    end
    % set up a frame to plot the iamge in
    p.stimWidthPix = p.stimHeightPix*size(im,2)/size(im,1);
    p.framePos=[p.centerPix(1)-p.stimWidthPix/2,p.centerPix(2)-p.stimHeightPix/2,p.centerPix(1)+p.stimWidthPix/2,p.centerPix(2)+p.stimHeightPix/2];
    
    % also load my prototype images which i'll show at start of each block.
    for ii=1:4
        
        imfn = fullfile(p.imagedir, sprintf('Shape_%.2f_%.2f.png', proto_coords(ii,1),proto_coords(ii,2)));
        if exist(imfn,'file')
            im=imread(imfn);
        else
            error('image file %s not found!',imfn)
        end        
        
        proto_ims(ii).name=imfn;
        proto_ims(ii).imtext=Screen('MakeTexture',w,im);
        p.proto_imfns{ii} = imfn;

    end
    
    % set up a set of frame to plot these images in
    p.protoStimWidthPix = p.protoStimHeightPix*size(im,2)/size(im,1);
   
    protocenters_pix = repmat(p.centerPix,4, 1) + p.protoDistPix*[1,-1; -1,-1; -1,1; 1,1];
    
    % the spatial position we present the images in is different every
    % block, but the labels always have the same correspondence.
    proto_orders = repmat(1:4, 8,1);
    proto_orders(5:8,:) = fliplr(proto_orders(5:8,:));
    for pp=2:4
       proto_orders([pp,pp+4],:) = circshift(proto_orders([pp,pp+4],:),pp-1,2);
    end
    
    protocenters_pix = protocenters_pix(proto_orders(datasample(1:8,1),:),:);
    p.protoFramePos = [ protocenters_pix(:,1)-p.protoStimWidthPix/2, ...
                        protocenters_pix(:,2)-p.protoStimHeightPix/2,...
                        protocenters_pix(:,1)+p.protoStimWidthPix/2,...
                        protocenters_pix(:,2)+p.protoStimHeightPix/2];
    p.protoTextPos = [ protocenters_pix(:,1), (protocenters_pix(:,2)+(p.protoStimHeightPix/2 + p.protoStimHeightPix/5))];
   %% keys
    KbName('UnifyKeyNames')

    %use number pad - change this for scanner 
    if scannerlaptop
        p.keys=[KbName('b'),KbName('y'),KbName('g'),KbName('r')];
    else
        p.keys =[KbName('1!'),KbName('2@'),KbName('3#'),KbName('4$')];
    end
    
    p.escape = KbName('escape');
    p.space = KbName('space');
    p.start = KbName('t');
    
      
   %% gamma correction
   gammacorrect = false;
   OriginalCLUT = [];
   
%     if gammacorrect
%         OriginalCLUT = Screen('LoadClut', w);
%         MyCLUT = zeros(256,3); MinLum = 0; MaxLum = 1;
%         CalibrationFile = 'calib_07-Nov-2016.mat';
%         [gamInverse,dacsize] = LoadCalibrationFileRR(CalibrationFile, expdir, GeneralUseScripts);
%         LumSteps = linspace(MinLum, MaxLum, 256)';
%         MyCLUT(:,:) = repmat(LumSteps, [1 3]);
%         MyCLUT = round(map2map(MyCLUT, repmat(gamInverse(:,4),[1 3]))); %Now the screen output luminance per pixel is linear!
%         Screen('LoadCLUT', w, MyCLUT);
%         clear CalibrationFile gamInverse
%     end
   
    %% Allocate arrays to store trial info
    data.Response = nan(p.nTrials,1); % response 1-4
    t.RespTimeFromOnset = nan(p.nTrials,1); % response time 
    
    %% timing information
    
    t.StimTimeTotal= 1.0;     
    % how long do they have to respond from image onset?
    t.RespTimeRange = 2;
    t.RespTimeBlank = t.RespTimeRange - t.StimTimeTotal;
    t.ShowFeedbackTime = 1;
    if p.Training
        t.MaxRespTime = 10;
    end
    
    t.ITIrange = [2,2];
    itis = linspace(t.ITIrange(1),t.ITIrange(2),p.nTrials);
    t.ITI = itis(randperm(length(itis)))';
   
    if scannerlaptop
        t.StartFixation = 13;
        t.EndFixation = 8;
    else
        t.StartFixation = 2;
        t.EndFixation = 0;
    end
    
    t.stim_flips = nan(p.nTrials,2);  
     %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    Screen(w,'TextFont','Helvetica');

    
%     InstrText = ['Quadrant Report Task.' '\n\n'...
%                 'Press buttons 1-4 to indicate quadrant.'];
    
    % draw the quadrant prototype images
    for qq=1:4
        Screen('DrawTexture', w, proto_ims(qq).imtext,[],p.protoFramePos(qq,:));
        DrawFormattedText(w, sprintf('%d',qq), p.protoTextPos(qq,1),p.protoTextPos(qq,2),p.white);
    end
            
%     DrawFormattedText(w, InstrText, 'center', 'center', p.white);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
    Screen('DrawingFinished', w);
    Screen('Flip', w);

    resp=0;
    % wait for a space bar press to start
    while resp==0
        [resp, timeStamp] = checkForResp([p.start,p.space],p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
%         [keyIsDown, secs, keyCode] = KbCheck([-1]);       
    end
    t.StartTime = GetSecs;
    
    KbReleaseWait();
    
    %% Fixation period before starting the stimuli (for scanner, this is the 12.8 seconds thing)
    FlushEvents('keyDown');
    Screen('Flip', w);
    ListenChar(2)
    
%     t.StartTime = GetSecs; %Get the starttime of the experiment in seconds
    GlobalTimer = 0; %this timer keeps track of all the timing in the experiment. TOTAL timing.
    TimeUpdate = t.StartTime; %what time is it now?
    % presentt begin fixation
%     Screen('BlendFunction', w, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    %TIMING!:
    GlobalTimer = GlobalTimer + t.StartFixation;
    TimePassed = (GetSecs-TimeUpdate); %Flush the time the previous event took
    while (TimePassed<t.StartFixation) %For as long as the cues are on the screen...
        TimePassed = (GetSecs-TimeUpdate);%And determine exactly how much time has passed since the start of the expt.       
        [resp, ~] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end; 
    end
    TimeUpdate = TimeUpdate + t.StartFixation;

    %% start trial loop

    for tt=1:p.nTrials

        %% Show target image   
        % start checking responses as soon as stim comes up
        keepChecking = 1;
        
        Screen('DrawTexture', w, allims(p.imlist(tt,1)).imtext,[],p.framePos);
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);

        t.stim_flips(tt,1) = GetSecs;

        Screen('Flip', w); 

        GlobalTimer = GlobalTimer + t.StimTimeTotal;

        TimePassed = (GetSecs-TimeUpdate);
        while TimePassed < t.StimTimeTotal
            TimePassed = (GetSecs-TimeUpdate);
            %check for escape responses 
            [resp, timeStamp] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end;
            if keepChecking && resp && find(p.keys==resp) 
                %they responded to this stim with 1-4
                data.Response(tt)=find(p.keys==resp);
                t.RespTimeFromOnset(tt)= timeStamp - t.stim_flips(tt,1);
                % now we have one response - stop checking! they can't
                % change it now even if they want to
                keepChecking=0;
            end
        end

        TimeUpdate = TimeUpdate + t.StimTimeTotal; %Update Matlab on what time it is.
       
        %% extra response screen to show prototypes if training block.
        if p.Training
            % draw the quadrant prototype images - in a random spatial
            % order but with correct labels.
            randseq = proto_orders(datasample(1:8,1),:);
            for qq=1:length(randseq)
                Screen('DrawTexture', w, proto_ims(qq).imtext,[],p.protoFramePos(randseq(qq),:));
                DrawFormattedText(w, sprintf('%d',qq), p.protoTextPos(randseq(qq),1),p.protoTextPos(randseq(qq),2),p.white);
            end
            
%             Screen('DrawTexture', w, allims(p.imlist(tt,1)).imtext,[],p.framePos);
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);

            t.stim_flips(tt,2) = GetSecs;

            Screen('Flip', w);
            
            TimePassed = (GetSecs-TimeUpdate);
            while keepChecking && TimePassed<t.MaxRespTime
                TimePassed = (GetSecs-TimeUpdate);
                %check for escape responses 
                [resp, timeStamp] = checkForResp(p.keys,p.escape);
                if resp==-1; escaperesponse(OriginalCLUT); end;
                if keepChecking && resp && find(p.keys==resp) 
                    %they responded to this stim with 1-4
                    data.Response(tt)=find(p.keys==resp);
                    t.RespTimeFromOnset(tt)= timeStamp - t.stim_flips(tt,2);
                    % now we have one response - stop checking! they can't
                    % change it now even if they want to
                    keepChecking=0;
                end
            end
            if keepChecking
                keepChecking = 0;
                data.Response(tt) = 0;
            end
            GlobalTimer = GlobalTimer + TimePassed;
            TimeUpdate = TimeUpdate + TimePassed; %Update Matlab on what time it is.
            
            %% Feedback (training only)
                       
            for qq=1:length(randseq)
                Screen('DrawTexture', w, proto_ims(qq).imtext,[],p.protoFramePos(randseq(qq),:));
                DrawFormattedText(w, sprintf('%d',qq), p.protoTextPos(randseq(qq),1),p.protoTextPos(randseq(qq),2),p.white);               
            end
            
             
            if data.Response(tt)~=0
                % outline their answer too
                Screen('FrameRect', w, p.fixColor, p.protoFramePos(randseq(data.Response(tt)),:),p.rectWidthPix);
            
                if data.Response(tt)==p.quadrant(tt)
                    showtext = 'Correct!\n\n';
                else
                    showtext = 'Incorrect.\n\nCorrect response is shown in green.';
                end
            else
               showtext = 'Time out!\n\nCorrect response is shown in green'; 
            end
            
            % outline the correct answer
            Screen('FrameRect', w, p.rectColor, p.protoFramePos(randseq(p.quadrant(tt)),:),p.rectWidthPix);
              
            
            DrawFormattedText(w, showtext, 'center', 'center', p.white);
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);
            Screen('Flip', w);
            
            GlobalTimer = GlobalTimer + t.ShowFeedbackTime;

            
            TimePassed = (GetSecs-TimeUpdate);
            while TimePassed<t.ShowFeedbackTime
                TimePassed = (GetSecs-TimeUpdate);
                %check for escape responses 
                [resp, timeStamp] = checkForResp(p.keys,p.escape);
                if resp==-1; escaperesponse(OriginalCLUT); end;
            end
            
            TimeUpdate = TimeUpdate + t.ShowFeedbackTime; %Update Matlab on what time it is.
                    
        end
        
        %% back to blank screen - keep checking for 1 second

        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);                     
        Screen('Flip', w);
        if ~p.Training
            t.stim_flips(tt,2)= GetSecs;
        end
        
        %TIMING!:
        GlobalTimer = GlobalTimer + t.RespTimeBlank;
        TimePassed = (GetSecs-TimeUpdate); 
        while TimePassed < t.RespTimeBlank
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, timeStamp] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end;            
            if keepChecking && resp && find(p.keys==resp)
                if keepChecking && resp && find(p.keys==resp)
                    %they responded to this stim with 1-4
                    data.Response(tt)=find(p.keys==resp);
                    t.RespTimeFromOnset(tt)= timeStamp - t.stim_flips(tt,1);
                    % now we have one response - stop checking! they can't
                    % change it now even if they want to
                    keepChecking=0;
                end
            end
            
           
            
        end
        % if they still haven't responded - this is when we would mark
        % it as a missed response.
        if keepChecking 
             keepChecking=0;
             data.Response(tt) = 0;
        end
        TimeUpdate = TimeUpdate + t.RespTimeBlank; %Update Matlab on what time it is.
           
        %% Now feedback
        
        if ~p.Training
            if data.Response(tt)~=0
                if data.Response(tt)==p.quadrant(tt)
                    showtext = 'Correct!\n\n';
                else
                    showtext = sprintf('Incorrect.\n\nCorrect answer is %d.',p.quadrant(tt));
                end
            else
               showtext = sprintf('Time out!\n\nCorrect answer is %d.',p.quadrant(tt));
            end
            
            DrawFormattedText(w, showtext, 'center', 'center', p.white);
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);
            Screen('Flip', w);
            
            GlobalTimer = GlobalTimer + t.ShowFeedbackTime;
            
            TimePassed = (GetSecs-TimeUpdate);
            while TimePassed<t.ShowFeedbackTime
                TimePassed = (GetSecs-TimeUpdate);
                %check for escape responses 
                [resp, timeStamp] = checkForResp(p.keys,p.escape);
                if resp==-1; escaperesponse(OriginalCLUT); end;
            end
            
            TimeUpdate = TimeUpdate + t.ShowFeedbackTime; %Update Matlab on what time it is.
                    
        end
        
        %% ITI 

        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);                     
        Screen('Flip', w);
       
        %TIMING!:
        GlobalTimer = GlobalTimer + t.ITI(tt,1);
        TimePassed = (GetSecs-TimeUpdate); 
        while TimePassed < t.ITI(tt,1)
            TimePassed = (GetSecs-TimeUpdate); 
            [resp, timeStamp] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end;            
        end
        
        TimeUpdate = TimeUpdate + t.ITI(tt,1); %Update Matlab on what time it is.
           
    end 
    %% finish experiment 
       
    % final fixation:
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    GlobalTimer = GlobalTimer + t.EndFixation;
    TimePassed = GetSecs-TimeUpdate;
    while (TimePassed<t.EndFixation) 
         TimePassed = GetSecs-TimeUpdate; 
        [resp, ~] = checkForResp(p.keys,p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end; 
    end
    
    t.EndTime = GetSecs; %Get endtime of the experiment in seconds
    t.TotalExpTime = (t.EndTime-t.StartTime); %Gets the duration of the total run.
    t.TotalExpTimeMins = t.TotalExpTime/60; %TOTAL exp time in mins including begin and end fixation.
    t.GlobalTimer = GlobalTimer; %Spits out the exp time in secs excluding begin and end fixation.

    %% get accuracy
    
    % don't worry about time out trials or the first trial
    inds2use = ~isnan(data.Response);
    
    acc = mean(data.Response(inds2use)==p.quadrant(inds2use));
   
    data.Accuracy = acc;
    
    fprintf('\nCompleted block %d!\n',p.runNumGlobal);
    fprintf('Percent accuracy: %.2f\n', acc*100);
    fprintf('Number of time out trials: %d/%d\n',sum(data.Response==0),p.nTrials);
    
    InstrText = ['Block finished!' '\n\n'...
                'Overall Accuracy: ' sprintf('%.2f',acc*100) '\n\n'];

    DrawFormattedText(w, InstrText, 'center', 'center', p.white);
    % put up a message to wait
    Screen('DrawingFinished', w);
    Screen('Flip', w);
         
    %----------------------------------------------------------------------
    %SAVE OUT THE DATA-----------------------------------------------------
    %----------------------------------------------------------------------
    if exist(p.fnsave_local,'file')
       load(p.fnsave_local);
    end
    
    %First I make a list of variables to save:
    TheData(p.runNumGlobal).t = t;
    TheData(p.runNumGlobal).p = p;
    TheData(p.runNumGlobal).data = data;
 
    save(p.fnsave_local,'TheData');
    if saveremote
        save(p.fnsave_remote,'TheData');
    end
    resp=0;
    % wait for a space bar press to exit
    while resp~=p.space
        [resp, timeStamp] = checkForResp(p.space, p.escape);
        if resp==-1; escaperesponse(OriginalCLUT); end 
    end

    KbReleaseWait();
    
    
    
    %----------------------------------------------------------------------
    %WINDOW CLEANUP--------------------------------------------------------
    %----------------------------------------------------------------------
    %This closes all visible and invisible screens and puts the mouse cursor
    %back on the screen
    Screen('CloseAll');
    if exist('OriginalCLUT','var')
        if exist('ScreenNumber','var')
            Screen('LoadCLUT', ScreenNumber, OriginalCLUT);
        else
            Screen('LoadCLUT', 0, OriginalCLUT);
        end
    end
    clear screen
    ListenChar(1);
    ShowCursor;
    
catch err%If an error occurred in the "try" block, this code is executed
   
    if exist('OriginalCLUT','var') && ~isempty(OriginalCLUT)
        if exist('ScreenNumber','var')
            Screen('LoadCLUT', ScreenNumber, OriginalCLUT);
        else
            Screen('LoadCLUT', 0, OriginalCLUT);
        end
    end
    Screen('CloseAll');                
    ShowCursor;
    if IsWin
        ShowHideWinTaskbarMex;     
    end
    ListenChar(1)
    rethrow(err)
%     if exist('ThrowErrorDB','file') ~= 0 %If ThrowErrorDB exists, use it
%         ThrowErrorDB; %Display last error (in a pretty way)
%     else
% %          rethrow(err)
%         disp('An error occured, but ThrowErrorDB is not in path, so the error cannot be displayed.');
%     end
end

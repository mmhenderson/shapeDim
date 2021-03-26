 % QUADRANT REPORT TASK - SHAPE SPACE 
% Subject sees a sequence of images - must report the categorical
% "quadrant" of each presented image. 
% Each quadrant has a corresponding prototype which is used to remind the
% subject of the relevant quadrants. 
 
% WRITTEN FOR BEHAVIOR ROOMS (B,C,D)

% parameters
    % Image size 24 degrees (fills almost entire screen)
    % Image onscreen for 1 second
    % Response period for max 10 second
    % Feedback 1 seconds
    % ITI 2 seconds
    % 48 trials if training==0, 24 trials if training==1  
try
    
    echo off
    clear 
    close all hidden
       
    p.scannerlaptop = 0;
   
    expdir = pwd;
    filesepinds = find(expdir==filesep);
    root = expdir(1:filesepinds(end)-1);
    
    % set up paths for saving
    if p.scannerlaptop || IsWindows
        % if scanner laptop or debugging, save to a subfolder here.
        datadir_local = fullfile(root, 'Data');
        datadir_remote = datadir_local;
    else
        % linux, behavior rooms - save to a different folder.
        datadir_local = '/home/pclexp/Documents/Maggie/shapeDim/Data/';
        datadir_remote = '/mnt/pclexp/Maggie/shapeDim/Data/';
    end
    
    %% Collect information about the subject, the date, etc.
    
    % make a dialog box to enter things   
    prompt = {'Debug mode?','Subject Initials','Subject Number','Session (1 or 2)',...
        'Run Number (1+)', 'Training Run? (0 or 1)','Image Set'};
    dlgtitle = 'Enter Run Parameters';
    dims = [1 35];
    definput = {'0','XX','99','1','1','0','3'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    p.debug = str2double(answer{1});
    p.SubInit = answer{2};
    p.SubNum = str2double(answer{3});
    p.SubNumStr = sprintf('%02d',str2double(answer{3}));
    p.Session = str2double(answer{4});
    p.RunNum = str2double(answer{5});
    p.Training = str2double(answer{6});
    p.ImageSet = str2double(answer{7});

    % check values of these entries
    if ~ismember(p.Session,[1,2]);  error('session must be 1 or 2');  end  
    if ~ismember(p.Training,[0,1]); error('training must be 0 or 1'); end
    if ~ismember(p.ImageSet,1:4); error('image set must be 1-4'); end

    rng('default')
    t.MySeed = sum(100*clock); 
    rng(t.MySeed);
    
    % where will i look for images?
%     p.imagedir=fullfile(root, sprintf('Stimuli/AmpGrid%d_full_grey/',p.ImageSet));
    p.imagedir=fullfile(root, sprintf('Stimuli/AmpGrid%d_adj_full_grey/',p.ImageSet));
      
    % am i debugging my code (usually on the admin account?) if so, don't
    % try to save remotely yet
    if p.debug
        saveremote = 0;
        Screen('Preference','SkipSyncTests',1);
    else
        saveremote = 1;
    end
    saveremote = 0; % never save remotely - running on lab computer w/o network connection 3/19/2021
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
         
    if p.Training
        p.expName = sprintf('PreTask_sess%d_TRAINING',p.Session);
    else
        p.expName = sprintf('PreTask_sess%d',p.Session);
    end
    
    p.fnsave_local = fullfile(datadir_local, ['S', p.SubNumStr, '_' p.expName '_' t.TheDate '.mat']);
    p.fnsave_remote = fullfile(datadir_remote, ['S', p.SubNumStr, '_' p.expName '_' t.TheDate '.mat']);
    
    if exist(p.fnsave_local,'file')
        load(p.fnsave_local);
        if p.RunNum~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.RunNum~=1            
        error('No data exists yet for this subject, check your run number')
    end
    
    p.rndseed = round(sum(100*clock));
    rng(p.rndseed);
    
    if p.debug
        % debug mode - do a few trials then quit
        p.nTrials = 8 ;
    elseif p.Training
        p.nTrials = 24;
    else
        p.nTrials = 48;
    end
 
    p.nIms = p.nTrials;
   
    %% figure out what to do during this block.
    
    % define a fixed list of all the possible finger-shape mappings that we
    % can use. We'll cycle through these based on subject number, keeping
    % the mapping the same for all runs of any given subject. 
    all_mappings = fliplr(perms(1:4));
    my_ind = mod(p.SubNum, size(all_mappings, 1));
    if my_ind==0
        my_ind = size(all_mappings,1);
    end
    p.mapping = all_mappings(my_ind,:);
%     p.mapping =[1,2,3,4];
    % this will also determine the order the shapes are shown onscreen
    % during the response screen (training only).
    
   %% set up the image grid 
   
   % first I'm defining all possible images that we can use in this task.
    start = 0;    % min value along each axis
    stop = 5; % max value along each axis  
    step = 0.1;
    center = (stop-start)./2+start;
      
    all_pts = round(start:step:stop,1);  
    [gridx,gridy] = meshgrid(all_pts,all_pts);
    all_grid_points = [gridx(:),gridy(:)];
    
    % now taking out images at exactly the prototype locations so that we
    % never use these during task
    nsteps_main = 4;
    start_grid = 0.1;
    stop_grid = 4.9;
    main_pts = round(linspace(start_grid,stop_grid, nsteps_main),1);
    
    proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
    proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
    proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
    assert(numel(proto_inds)==4)
    all_grid_points(proto_inds,:) = [];
  
    % Also taking out any images along quadrant boundaries because these
    % are ambiguous 
    bound_inds = find(any(all_grid_points==center,2));
    all_grid_points(bound_inds,:) = [];
  
    % now define which quadrant each image lies in (this is a physical
    % property that is the same regardless of what task is being done).
    % Category (defined below) will change with task.
    all_quadrant = zeros(size(all_grid_points,1),1);
    all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)>center) = 1;
    all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)>center) = 2;
    all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)<center) = 3;
    all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)<center) = 4;
    
    % Next, for each point in the full grid, define the difficulty level 
    % based on distance from the boundary
  
    % considering distance from any boundary
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

    nperbin = sum(repmat(dist_groups_binned, 1, size(bin_edges,2))==repmat(1:size(bin_edges,2),size(dist_groups_binned,1),1),1);
    assert(all(nperbin>=p.nIms/4));
    
    %% Now choose the images that we want to show on this block. Sampling from multiple difficulty levels.
    if p.Training
        difficulty_to_sample = 6:9;
    else
        difficulty_to_sample = 6:13;
    end
    p.difficulty_to_sample = difficulty_to_sample;
    assert(~mod(p.nIms, length(difficulty_to_sample)));
    nPerDiff = p.nIms/length(difficulty_to_sample);
    myinds = [];
    for dd=1:length(difficulty_to_sample)
        thisbin = find(dist_groups_binned==difficulty_to_sample(dd));    
        myinds = [myinds; datasample(thisbin, nPerDiff, 'replace',false)];        
    end
  
    p.points = all_grid_points(myinds,:);
    p.difficulty = dist_from_bound(myinds);
    p.quadrant = all_quadrant(myinds);
     
    % check this math
    dist_from_bound = round(min(abs((p.points-repmat(center, size(p.points,1),2))),[],2),1);
    assert(all(dist_from_bound==repelem(bin_edges(difficulty_to_sample)',nPerDiff))); 
        
    %% Make a randomized image sequence    
    p.imlist = (1:p.nIms)';
    p.imlist = p.imlist(randperm(length(p.imlist)));
    
    % sort everything the same way
    p.points = p.points(p.imlist,:);
    p.difficulty= p.difficulty(p.imlist);
    p.quadrant = p.quadrant(p.imlist);
    
    %% define the category for each image
    % p.category changes with the mapping - p.quadrant is constant no matter the
    % task being performed or the mapping.
    p.category = zeros(size(p.quadrant));
    for ii=1:4
        p.category(p.quadrant==ii) = p.mapping(ii);
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
%     plot(p.points(p.category==1,1),p.points(p.category==1,2),'ro');
%     plot(p.points(p.category==2,1),p.points(p.category==2,2),'bo');
%     plot(p.points(p.category==3,1),p.points(p.category==3,2),'go');
%     plot(p.points(p.category==4,1),p.points(p.category==4,2),'ko');
%     line([center,center],get(gca,'YLim'),'Color','k')
%     line(get(gca,'XLim'), [center,center], 'Color','k');
%     
%     set(gcf,'Color','w');
%     title(sprintf('Shape locations used for this block: Quadrant Task\n'));
%     
   
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
        p.backColor = 77;   % round(0.3*255)
        if p.debug
            Screen('Preference', 'SkipSyncTests', 1);
        end
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

    if p.scannerlaptop
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
        p.vDistCM = 46; 
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
    p.stimHeightDeg = 24;   
    % note that this is the size of the square image that is drawn, including its grey background. 
    % The shape itself is about 2/3 of that size. So, the background is
    % drawn in a frame slightly bigger than the screen, but all the shape pixels are
    % within the bounds of the screen.

    p.protoStimHeightDeg = 9;   % how big is each prototype image?
%     p.protoSpacingDeg = 1;  % how far apart are the four prototypes drawn?
     
    p.outlineProp = 0.9; % for the colored square outlines, how big are they relative to the image square size?
    p.feedbackTextHeightDeg = 8;
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
    % set up a frame to plot the image in
    p.stimWidthPix = p.stimHeightPix*size(im,2)/size(im,1);
    p.framePos=[p.centerPix(1)-p.stimWidthPix/2,p.centerPix(2)-p.stimHeightPix/2,p.centerPix(1)+p.stimWidthPix/2,p.centerPix(2)+p.stimHeightPix/2];
    
    % also load my prototype images which i'll show at start of each block.
    % note that we are always loading the same images in the same order to
    % put inside the proto_ims structure, but we will choose the response
    % mapping differently depending on what the experimenter entered.
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
    
   
    p.protoStimWidthPix = p.protoStimHeightPix*size(im,2)/size(im,1);
    p.boundTextOutlinePix = p.protoStimWidthPix/6;
    %% Define screen positions for plotting the prototype images
    smalldist = p.protoStimWidthPix/2+p.protoSpacingPix;
    largedist = smalldist + p.protoSpacingPix + p.protoStimWidthPix;
    protocenters_pix = repmat(p.centerPix,4,1) + [-largedist,0;-smalldist,0;smalldist,0;largedist,0];
%     protocenters_pix = repmat(p.centerPix,4, 1) + p.protoDistPix*[-3,0; -1,0; 1,0; 3 ,0];    
   
    % define where to draw the images themselves - this is a frame to pass
    % into the drawtexture function
    p.protoFramePos = [ protocenters_pix(:,1)-p.protoStimWidthPix/2, ...
                        protocenters_pix(:,2)-p.protoStimHeightPix/2,...
                        protocenters_pix(:,1)+p.protoStimWidthPix/2,...
                        protocenters_pix(:,2)+p.protoStimHeightPix/2];
                    
    % now create two outline boxes, one for each category
    % first list all four outlines (to help draw the boxes)
    outlines = [ protocenters_pix(:,1)-p.outlineProp*p.protoStimWidthPix/2, ...
                        protocenters_pix(:,2)-p.outlineProp*p.protoStimHeightPix/2,...
                        protocenters_pix(:,1)+p.outlineProp*p.protoStimWidthPix/2,...
                        protocenters_pix(:,2)+p.outlineProp*p.protoStimHeightPix/2];
    p.protoOutlinePos = outlines;
  
    %% keys
    KbName('UnifyKeyNames')

    %use number pad - change this for scanner 
    if p.scannerlaptop
        p.keys=[KbName('b'),KbName('y'),KbName('g'),KbName('r')];
    else
%         p.keys =[KbName('1!'),KbName('2@'),KbName('3#'),KbName('4$')];
        p.keys=[KbName('u'),KbName('i'),KbName('o'),KbName('p')];
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
%     if p.Training
    t.MaxRespTime = 10;
%     end 
    
    t.ITIrange = [2,2];
    itis = linspace(t.ITIrange(1),t.ITIrange(2),p.nTrials);
    t.ITI = itis(randperm(length(itis)))';
   
%     if p.scannerlaptop
%         t.StartFixation = 13;
%         t.EndFixation = 8;
%     else
    t.StartFixation = 2;
    t.EndFixation = 0;
%     end
    
    t.stim_flips = nan(p.nTrials,2);  
     %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    if p.scannerlaptop==1
        Screen(w,'TextFont','Helvetica');
    else
        Screen(w,'TextFont','-bitstream-courier 10 pitch-medium-i-normal--0-0-0-0-m-0-adobe-standard')
    end
    
    % draw the prototype images
    for bb=1:4
        ind = find(p.mapping==bb);
        Screen('DrawTexture', w, proto_ims(ind).imtext,[],p.protoFramePos(bb,:));   
%         Screen('DrawTexture', w, proto_ims(bb).imtext,[],p.protoFramePos(p.mapping(bb),:));        
%         DrawFormattedText(w, sprintf('%d',bb), p.boundTextPos(bb,1)+p.textAdjPix(1),p.boundTextPos(bb,2)+p.textAdjPix(2),p.white);
    end
  
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
        
        Screen('DrawTexture', w, allims(tt).imtext,[],p.framePos);
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
       
        %% put up the prototypes to remind the subject of responses

        % draw the prototype images
        for bb=1:4
            Screen('DrawTexture', w, proto_ims(bb).imtext,[],p.protoFramePos(p.mapping(bb),:));        
%                 DrawFormattedText(w, sprintf('%d',bb), p.boundTextPos{oo}(bb,1)+p.textAdjPix(1),p.boundTextPos{oo}(bb,2)+p.textAdjPix(2),p.white);
        end

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
                t.RespTimeFromOnset(tt)= timeStamp - t.stim_flips(tt,1);
                % now we have one response - stop checking! they can't
                % change it now even if they want to
                keepChecking=0;
            end
        end
        % if they still haven't responded - this is when we would mark
        % it as a missed response.
        if keepChecking
            keepChecking = 0;
            data.Response(tt) = 0;
        end
        GlobalTimer = GlobalTimer + TimePassed;
        TimeUpdate = TimeUpdate + TimePassed; %Update Matlab on what time it is.

        %% Feedback 

       % draw the quadrant prototype images
        for bb=1:4
            ind = find(p.mapping==bb);
            Screen('DrawTexture', w, proto_ims(ind).imtext,[],p.protoFramePos(bb,:));        
%                 DrawFormattedText(w, sprintf('%d',bb), p.boundTextPos{oo}(bb,1)+p.textAdjPix(1),p.boundTextPos{oo}(bb,2)+p.textAdjPix(2),p.white);
        end

        if data.Response(tt)~=0
             % outline their answer
            Screen('FrameRect', w, p.fixColor, p.protoOutlinePos(data.Response(tt),:),p.rectWidthPix);

            if data.Response(tt)==p.category(tt)
                showtext = 'Correct!\n\n';
            else
                showtext = 'Incorrect.\n\nCorrect response is shown in green.';
            end
        else 
           showtext = 'Time out!\n\nCorrect response is shown in green.'; 
        end
%             showtext = sprintf('True quadrant was %d, true response was %d',p.quadrant(tt),p.category(tt));
        % outline the correct answer
        Screen('FrameRect', w, p.rectColor, p.protoOutlinePos(p.category(tt),:),p.rectWidthPix);

        DrawFormattedText(w, showtext, 'center', p.centerPix(2)-p.feedbackTextHeightPix, p.white);
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
       
    
    % get accuracy
    trialsdone = ~isnan(t.stim_flips(:,2));
    acc = mean(data.Response(trialsdone)==p.category(trialsdone));
    data.Accuracy = acc;
    
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
    
    fprintf('\nCompleted block %d!\n',p.RunNum);
    fprintf('Accuracy is %.2f percent\n',data.Accuracy*100);
    fprintf('Number of time out trials: %d/%d\n',sum(data.Response==0),p.nTrials);
    
    InstrText = ['Block finished!' '\n\n'...
                sprintf('Accuracy is %.2f percent',data.Accuracy*100), '\n\n'...
                'Please find the experimenter to start next block.'];

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
    TheData(p.RunNum).t = t;
    TheData(p.RunNum).p = p;
    TheData(p.RunNum).data = data;
 
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

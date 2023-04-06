            % ONE BACK REPEAT DETECTION - SHAPE SPACE 
% Subject sees a sequence of images - must report whenever they see a
% "repeat" image. Press 1 for repeat, 2 for different.
 
% WRITTEN FOR BEHAVIOR ROOMS AND SCANNER LAPTOP
% change variable p.scannerlaptop to indicate which situation.
       
% parameters for MRI task
    % Image size 24  degrees (fills almost entire screen)
    % Image onscreen for 1 second    
    % Response period for 1 second (blank)
    % ITI [1-5] seconds jittered
    
    % 48 total images presented
    
    % total time = 13s blank start + 48*(1+1+3) + 10s blank end 
    %   = 263 seconds or 4:23
    %   with TR=0.8 sec this is 329 TRs.
    
try
    
    echo off
    clear 
    close all hidden
       
    p.scannerlaptop = 1 ;
   
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
    prompt = {'Debug mode?','Subject Initials','Subject Number','Session (1-3)',...
       'Run Number','Image Set','Difficulty (1-7)','Low-contrast version?'};
    dlgtitle = 'Enter Run Parameters';
    dims = [1 35];
    definput = {'0','BX','2','3','1','3','3','1    '};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    p.debug = str2double(answer{1});
    p.SubInit = answer{2};  
    p.SubNum = str2double(answer{3 });
    p.SubNumStr = sprintf('%02d',str2double(answer{3}));
    p.Session = str2double(answer{4});
    p.RunNum = str2double(answer{5});
    p.ImageSet = str2double(answer{6});
    p.RunDifficulty = str2double(answer{7});
    p.LowContrastVersion = boolean(str2double(answer{8}));
    
    % check values of these entries
    if ~ismember(p.Session,[1,2,3]);  error('session must be 1-3');  end
    if ~ismember(p.ImageSet,[3]); error('image set must be 3'); end
    if ~ismember(p.RunDifficulty,1:7); error('difficulty level must be 1-7'); end
   
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
         
   
    p.expName = sprintf('OneBackTaskMRI_sess%d',p.Session);
    
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
    
    p.nTrials = 48;
    p.nTrials_variable = 16;    % variable are the trials whose difficulty is set by experimenter
    p.nTrials_main = 32;    % main are trials from the fixed 16-pt grid.

    if p.debug  
        % debug mode - do a few trials then quit
          p.nTrials = 8;      
    end

    % choose response mapping now. Always alternates even/odd runs. Whether
    % we start the session with mapping 1 or 2 will depend on subject
    % number (alternate odd/even).
    if (mod(p.RunNum,2)==1 && mod(p.SubNum,2)==1) ||...
            (mod(p.RunNum,2)==0 && mod(p.SubNum,2)==0)
        p.RespMap=1;
    else
        p.RespMap=2;
    end
    %% set up the main 16 pts grid 
   
    % fixed number of images in the grid
    p.nIms_main = 16;
    % Define a 4x4 grid w even spacing
    start = 0;    % min value along each axis
    stop = 5; % max value along each axis  
    nsteps_main = 4;    % how many steps along each axis?
    step = 0.1;
    start_grid = 0.1;
    stop_grid = 4.9;
    main_pts = round(linspace(start_grid,stop_grid, nsteps_main),1);
    
    [gridx,gridy] = meshgrid(main_pts,main_pts);
    main_grid_points = [gridx(:),gridy(:)];
    p.main_grid_points = main_grid_points;
    
    %% set up the points to use for the difficult trials
    variable_points = zeros(p.nTrials_variable,2);
    max_dist = diff(main_pts(1:2))/2;
    % for each difficulty level, define the range of distances that fall
    % into that difficulty level. Rows are difficulty levels, top is easiest.
    difficulty_distance_bins = flipud([step:step:(max_dist-step); 2*step:step:max_dist]'-0.05);
    dist_range = difficulty_distance_bins(p.RunDifficulty,:);
    nDiffLevels = size(difficulty_distance_bins,1);
    
    for pp = 1:p.nTrials_variable
        
        % for this point, define all the points that neighbor it
        pt = main_grid_points(pp,:);        
        [nx,ny] = meshgrid(pt(1)-max_dist+step:step:pt(1)+max_dist-step, pt(2)-max_dist+step:step:pt(2)+max_dist-step);
        nx = round(nx,1);
        ny = round(ny,1);
        all_neighbors = [nx(:),ny(:)];
        % now f ilter out only the ones within the desired "distance" from
        % the point of interest.
        in_bounds = nx(:)<=stop & ny(:)<=stop & nx(:)>=start & ny(:)>=start;
        rad = sqrt(sum((all_neighbors-repmat(pt, size(all_neighbors,1),1)).^2,2));
        
        neighbors = all_neighbors(in_bounds & rad>=dist_range(1) & rad<dist_range(2),:);
       
        variable_points(pp,:) = neighbors(datasample(1:size(neighbors,1),1),:); 
    end
    
    %% plot these if you want to see how the grid looks
    
%     figure;hold all;
%     
%     [x,y] = meshgrid(start:step:stop, start:step:stop);
%     allpts = [x(:),y(:)];
%     center = (stop-start)./2+start;
%      
%     plot(allpts(:,1),allpts(:,2),'.','MarkerEdgeColor',[0.8, 0.8, 0.8], 'MarkerFaceColor',[0.8, 0.8, 0.8]);
%     axis equal
%     plot(p.main_grid_points(:,1),p.main_grid_points(:,2),'bo');
%     plot(variable_points(:,1),variable_points(:,2),'ro');
%    
%     cols = parula(nDiffLevels);
%     for pp = 1:size(main_grid_points,1)
%         
%         % for this point, define all the points that neighbor it
%         pt = main_grid_points(pp,:);        
%         [nx,ny] = meshgrid(pt(1)-max_dist+step:step:pt(1)+max_dist-step, pt(2)-max_dist+step:step:pt(2)+max_dist-step);
%         nx = round(nx,1);
%         ny = round(ny,1);
%         all_neighbors = [nx(:),ny(:)];
%         
%         in_bounds = nx(:)<=stop & ny(:)<=stop & nx(:)>=start & ny(:)>=start;
% 
%         for dd  = 1:nDiffLevels
%             % now filter out only the ones within the desired "distance" from
%             % the point of interest.
%             dist_range = difficulty_distance_bins(dd,:);
%             neighbors = all_neighbors(in_bounds & rad>=dist_range(1) & rad<dist_range(2),:);
%             plot(neighbors(:,1),neighbors(:,2),'.','MarkerEdgeColor',cols(dd,:), 'MarkerFaceColor',cols(dd,:))        
%         end
%     end    
% 
%     line([center,center],get(gca,'YLim'),'Color','k')
%     line(get(gca,'XLim'), [center,center], 'Color','k');
%     
%     set(gcf,'Color','w');
%     title(sprintf('Shape locations used for this block: Difficulty %d\n',p.RunDifficulty));
%     
    %% create the full image sequence
    
    % first, finding a sequence of the main images (each 2x) where there
    % are no repeats to start with. brute force randomization until we find
    % a good sequence.   
    main_im_list = repmat((1:p.nIms_main)',2,1);
    any_repeats = 1;
    max_iters=1000;
    iter=0;
    while any_repeats && iter<max_iters
        iter=iter+1;
        main_im_list = main_im_list(randperm(length(main_im_list)));
        any_repeats = any(diff(main_im_list)==0);
    end
    
    % next, going to add in the repeat trials where we want them.
    % start building lists of the properties of images to be shown on each
    % trial (including only the main grid pts for now)
    points = p.main_grid_points(main_im_list,:);
    variable_im_list = nan(size(main_im_list));   
    is_repeat = zeros(size(main_im_list));
    
    % of the 16 variable trials, decide now which of those will be
    % identical repeats and which will be similar but nonidentical
    variable_repeat_freq = 0.5;
    ims2repeat = linspace(0,1,p.nTrials_variable)<variable_repeat_freq;
    ims2repeat = ims2repeat(randperm(length(ims2repeat)));
    % now insert them into the sequence
    for pp = 1:p.nTrials_variable
        % for this main grid image, choose at random one of its
        % presentations to be followed by a difficult trial
        ind = datasample(find(main_im_list==pp),1);
        % decide whether to use an identical repeat, or the nearby image
        % chosen above
        if ims2repeat(pp)==1
            newpt = p.main_grid_points(pp,:);
            is_repeat = [is_repeat(1:ind); 1; is_repeat(ind+1:end)];
        else
            newpt = variable_points(pp,:);
            is_repeat = [is_repeat(1:ind); 0; is_repeat(ind+1:end)];
        end
        % insert the desired image into all the growing lists
        main_im_list = [main_im_list(1:ind); nan; main_im_list(ind+1:end)];
        variable_im_list = [variable_im_list(1:ind); pp; variable_im_list(ind+1:end)];
        points = [points(1:ind,:); newpt; points(ind+1:end,:)];
        
    end
    
    % can use this later to decide which trials to use for analysis
    p.is_main_grid = ~isnan(main_im_list); 
    p.points = points;
    p.is_repeat = is_repeat==1;
    p.repeat_prob = mean(p.is_repeat);
    
    % determine what the correct finger response will be - this depends on
    % the response mapping in place (p.RespMap)
    if p.RespMap==1
        % 1 for repeat, 2 for no
        p.correct_resp = 2-is_repeat;
    else
        % 2 for repeat, 1 for no
        p.correct_resp = 1+is_repeat;
    end
    % first trial doesn't require a response
    p.correct_resp(1) = 0;
    %% calculate difficulty (e.g. distance) of each trial from previous
    % this is mostly to check calculations, also might be useful later in
    % analysis
    p.dist_from_previous = zeros(size(p.is_repeat));
    for tt = 2:size(p.dist_from_previous,1)
       p.dist_from_previous(tt) = sqrt(sum((p.points(tt,:) - p.points(tt-1,:)).^2));
    end
    
    % make sure all "hard" trials are in the expected distance range
    dist_range = difficulty_distance_bins(p.RunDifficulty,:);
    vals2check = p.dist_from_previous(~p.is_main_grid);
    assert(all((vals2check>=dist_range(1) & vals2check<dist_range(2)) | vals2check==0)); % can also be zero
    % and make sure repeat/nonrepeat labels are correct
    assert(all(p.dist_from_previous(p.is_repeat)==0));
    assert(all(p.dist_from_previous(~p.is_repeat & (1:length(p.dist_from_previous) )'>1)>0)); % first trial is zero    
  
    % adjust length if in debug mode
    if p.debug==1 
        p.is_main_grid = p.is_main_grid(1:p.nTrials);
        p.points = p.points(1:p.nTrials,:);
        p.is_repeat = p.is_repeat(1:p.nTrials);
        p.correct_resp = p.correct_resp(1:p.nTrials);
        p.dist_from_previous = p.dist_from_previous(1:p.nTrials);
    end
        
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
%     p.backColor = 77;   % round(0.3*255)
    % stuff about colors/grays
    p.backColorOrig = 77;
    p.shapeColorOrig = 230;
    if p.LowContrastVersion
        p.backColor = 50;
        p.fixColor = [0.05, 0.05, 0.05  ] * 255;
        p.shapeColor = 0.12 * 255;
    else
        p.backColor = p.backColorOrig;
        p.fixColor= [0.80, 0.80, 0.80]*255;
        p.shapeColor = p.shapeColorOrig;
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
        p.screenHeightCM = 16;
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
    
    
    
%     p.fixColor = [0,0,0];
    p.stimHeightDeg = 24;   
    % note that this is the size of the square image that is drawn, including its grey background. 
    % The shape itself is about 2/3 of that size. So, the background is
    % drawn in a frame slightly bigger than the screen, but all the shape pixels are
    % within the bounds of the screen.

%     p.protoStimHeightDeg = 9;   % how big is each prototype image?
%     p.protoSpacingDeg = 1;  % how far apart are the four prototypes drawn?
%     p.outlineProp = 0.9; % for the colored square outlines, how big are they relative to the image square size?
%     p.feedbackTextHeightDeg = 8;
%     p.rectColor = [0.1,0.8,0.1]*255;
%     p.rectWidthDeg = 0.1;
    % convert from degrees to pixel units
    p = deg2pix(p);
    p.fixSizePix = ceil(p.fixSizePix);
    %% Load the images

    for ii=1:p.nTrials
        
        imfn = fullfile(p.imagedir, sprintf('Shape_%.2f_%.2f.png', p.points(ii,1),p.points(ii,2)));
        if exist(imfn,'file')
            im=imread(imfn);
        else
            error('image file %s not found!',imfn)
        end        
         if p.LowContrastVersion
            im(im(:)==p.shapeColorOrig) = p.shapeColor;
            im(im(:)==p.backColorOrig) = p.backColor;
        end
        allims(ii).name=imfn;
        allims(ii).imtext=Screen('MakeTexture',w,im);
        p.imfns{ii} = imfn;

    end
    % set up a frame to plot the image in
    p.stimWidthPix = p.stimHeightPix*size(im,2)/size(im,1);
    p.framePos=[p.centerPix(1)-p.stimWidthPix/2,p.centerPix(2)-p.stimHeightPix/2,p.centerPix(1)+p.stimWidthPix/2,p.centerPix(2)+p.stimHeightPix/2];
  
    %% keys
    KbName('UnifyKeyNames')

    %use number pad - change this for scanner 
    if p.scannerlaptop
        p.keys=[KbName('b'),KbName('y')];
    else
%         p.keys =[KbName('1!'),KbName('2@'),KbName('3#'),KbName('4$')];
        p.keys=[KbName('u'),KbName('i')];
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

    if p.scannerlaptop
        % actual timing for scanner expt - long delay at start,
        % longer/jittered ITIs
        t.StartFixation = 13;
        t.EndFixation = 10; % making this run just slightly longer, so it isn't identical in len to main task.
        t.ITIrange = [1, 5];
    else
        % for behav room, no need to have these longer delays so whole expt
        % will be shorter.
        t.StartFixation = 2;
        t.EndFixation = 0;
        t.ITIrange = [2, 2];
    end
    
    itis = linspace(t.ITIrange(1),t.ITIrange(2),p.nTrials);
    t.ITI = itis(randperm(length(itis)))';
   
    t.stim_flips = nan(p.nTrials,2);  
     %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    if p.scannerlaptop==1
        Screen(w,'TextFont','Helvetica');
    else
        Screen(w,'TextFont','-bitstream-courier 10 pitch-medium-i-normal--0-0-0-0-m-0-adobe-standard')
    end

    if p.RespMap==1
        InstrText = char(['Repeat detection task.' '\n\n'...
                'Index finger = repeat.' '\n\n'...
                'Middle finger = no repeat.']);
    else
        InstrText = char(['Repeat detection task.' '\n\n'...
                'Middle finger = repeat.' '\n\n'...
                'Index finger = no repeat.']);
    end
    
    DrawFormattedText(w, InstrText, 'center', 'center', p.fixColor);
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
 
        %% Back to blank screen, but keep checking responses

        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);                     
        Screen('Flip', w);

        t.stim_flips(tt,2)= GetSecs;

        %TIMING!:
        GlobalTimer = GlobalTimer + t.RespTimeBlank;
        TimePassed = (GetSecs-TimeUpdate); 
        while TimePassed < t.RespTimeBlank
            TimePassed = (GetSecs-TimeUpdate); 
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
             keepChecking=0;
             data.Response(tt) = 0;
        end
        TimeUpdate = TimeUpdate + t.RespTimeBlank; %Update Matlab on what time it is.
      
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
    trialsdone(1) = 0;  % first trial didn't require a response, so don't count it here in accuracy calc
    acc = mean(data.Response(trialsdone)==p.correct_resp(trialsdone));
    dprime = get_dprime(data.Response(trialsdone), p.correct_resp(trialsdone));
    data.Dprime = dprime;
    data.Accuracy = acc; 
    
    data.MainGridAccuracy = mean(data.Response(trialsdone & p.is_main_grid==1)==...
        p.correct_resp(trialsdone & p.is_main_grid==1));
%     data.MainGridDprime = get_dprime(data.Response(trialsdone & p.is_main_grid==1),...
%         p.correct_resp(trialsdone & p.is_main_grid==1));
    data.VariableAccuracy = mean(data.Response(trialsdone & p.is_main_grid==0)==...
        p.correct_resp(trialsdone & p.is_main_grid==0));
    data.VariableDprime = get_dprime(data.Response(trialsdone & p.is_main_grid==0),...
        p.correct_resp(trialsdone & p.is_main_grid==0));
 
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
    fprintf('Accuracy is %.2f percent (dprime = %.2f)\n',data.Accuracy*100, data.Dprime);
    fprintf('Easier trials: %.2f percent\n',data.MainGridAccuracy*100);
    fprintf('Harder trials (diff=%d): %.2f percent (dprime = %.2f)\n',p.RunDifficulty,data.VariableAccuracy*100, data.VariableDprime);
    fprintf('Number of time out trials: %d/%d\n',sum(data.Response==0),p.nTrials);
    
    InstrText = ['Block finished!' '\n\n'...
                sprintf('Accuracy is %.2f percent (d-prime = %.2f)',data.Accuracy*100, data.Dprime)];

    DrawFormattedText(w, InstrText, 'center', 'center', p.fixColor);
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


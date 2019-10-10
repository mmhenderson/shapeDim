% ONE-BACK TASK - SHAPE SPACE 
% Subject sees a sequence of images presented one at a time at fixation.
% They perform repeat detection task, where they must respond whether the
% current image matches the previous one in sequence. 

% WRITTEN FOR SCANNER LAPTOP

% parameters
    % Image size 10 degrees (fills almost entire screen)
    % Image onscreen for 1 second
    % ITI jittered 2-8 seconds (5 average)
    
    % 48 total images presented (16 exemplars x 3 pres each)
    
    % total time = 13s blank start + 48*(1+5) + 8s end fixation = 309 s
    % 309 s/60 = 5.15 mins = 5:09
    
    % TR = 0.80 seconds
    % 309/0.8 = 387 TRs
    
try
%     
    echo off
    clear 
    close all hidden
       
    scannerlaptop = 1;
%     saveremote = 0;
    
    expdir = pwd;
    filesepinds = find(expdir==filesep);
    root = expdir(1:filesepinds(end)-1);
    datadir_local = fullfile(root, 'Data');

    %% Collect information about the subject, the date, etc.
    p.Subject = input('Initials? (default is tmp)--> ','s'); 
    if isempty(p.Subject)
        p.Subject = 'tmp'; 
    end 
    assert(strcmp(p.Subject,'CG'));
    
    SubNum = input('Subject number?  --> ','s'); 
    if isempty(SubNum)
        error('enter a valid subject number')
    end
    p.SubNum = sprintf('%02d', str2double(SubNum));
      
    RunNum = input('Run?  --> ','s');
    if isempty(RunNum)
        error('enter a valid run number')
    end
    p.runNumGlobal = str2double(RunNum);
         
    rng('default')
    p.rndseed = sum(100*clock); 
    rng(p.rndseed);

    p.expName = 'OneBackPilotTask';

    % where will i look for images? even runs and odd runs will switch.
    if mod(p.runNumGlobal,2)
        p.stimset = 3;
        p.imagedir=fullfile(root, 'Stimuli/AmpGrid3_grey/');
    else
        p.stimset = 4;
        p.imagedir=fullfile(root, 'Stimuli/AmpGrid4_grey/');
    end
    
    %% initialize my data file

    % save a copy of the currently running script here (in text format) so
    % we can look at it later if issues arise
    p.MyScript = fileread([mfilename('fullpath'),'.m']);
    
    % make sure my data dir exists
    if ~exist(datadir_local,'dir')
        mkdir(datadir_local)
    end
%     if saveremote && ~exist(datadir_remote,'dir')
%         mkdir(datadir_remote)
%     end
    t.TheDate = datestr(now,'yymmdd'); %Collect todays date (in t.)
    t.TimeStamp = datestr(now,'HHMM'); %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
    
    p.fnsave_local = fullfile(datadir_local, ['S', p.SubNum, '_' p.expName '_' t.TheDate '.mat']);
%     p.fnsave_remote = fullfile(datadir_remote, ['S', p.SubNum, '_' p.expName '_' t.TheDate '.mat']);
    
    if exist(p.fnsave_local,'file')
        load(p.fnsave_local);
        if p.runNumGlobal~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.runNumGlobal~=1            
        error('No data exists yet for this subject, check your run number')
    end
    
    p.nTrials = 48;
    p.nIms = 16;
    assert(mod(p.nTrials, p.nIms)==0)
    %% set up the grid 
   
    % Define a 4x4 grid w even spacing
    start = 0.2;    % min value along each axis
    stop = 4.8; % max value along each axis  
    nsteps_main = 4;    % how many steps along each axis?
    all_pts = round(linspace(start,stop, nsteps_main),1);
    
    [gridx,gridy] = meshgrid(all_pts,all_pts);
    all_grid_points = [gridx(:),gridy(:)];
    p.points = all_grid_points;
    %% Make a randomized image sequence
    
    imlist_unshuffled = repelem(1:p.nIms, p.nTrials/p.nIms)';
    
    p.imlist = imlist_unshuffled(randperm(length(imlist_unshuffled)));
    p.imcoords = all_grid_points(p.imlist,:);
    p.isrepeat = [0; diff(p.imlist)==0];
    p.correct_resp=2-p.isrepeat;
    p.num_repeats = sum(p.isrepeat);    
        
    %% set up my screen 
    
    InitializeMatlabOpenGL;  
    PsychImaging('PrepareConfiguration');
    Screen('Preference', 'SkipSyncTests', 0);
    AssertOpenGL; % bail if current version of PTB does not use
    PsychJavaTrouble;
        
    % don't change this, matches the images!
    p.backColor = 77;   % round(0.3*255)
    
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
%     p.stimHeightDeg = 14;
    p.stimHeightDeg = 24;   
    % note that this is the size of the square image that is drawn, including its grey background. 
    % The shape itself is about 2/3 of that size. So, the background is
    % drawn in a frame slightly bigger than the screen, but all the shape pixels are
    % within the bounds of the screen.
    
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

    %% keys
    KbName('UnifyKeyNames')

    %use number pad - change this for scanner 
    if scannerlaptop
%         p.keys=[KbName('b'),KbName('y'),KbName('g'),KbName('r')];
        p.keys=[KbName('b'),KbName('y')];
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
    
    t.ITIrange = [2,8];
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
    
    InstrText = sprintf('Fixate.\nPress 1 (index finger) for repeat.\nPress 2 (middle finger) for no repeat.\n\n\n\n');
    DrawFormattedText(w, InstrText, 'center', 'center', p.fixColor);
  
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
                %they responded to this stim with 1-2
                data.Response(tt)=find(p.keys==resp);
                t.RespTimeFromOnset(tt)= timeStamp - t.stim_flips(tt,1);
                % now we have one response - stop checking! they can't
                % change it now even if they want to
                keepChecking=0;
            end
        end

        TimeUpdate = TimeUpdate + t.StimTimeTotal; %Update Matlab on what time it is.
       
        %% back to blank screen for ITI - keep checking for 1 second

        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);                     
        Screen('Flip', w);
              
        t.stim_flips(tt,2) = GetSecs;
        
        %TIMING!:
        GlobalTimer = GlobalTimer + t.ITI(tt);
        TimePassed = (GetSecs-TimeUpdate); 
        while TimePassed < t.ITI(tt)
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
            % check if we're past the response range
            if keepChecking && TimePassed>(t.RespTimeRange-t.StimTimeTotal)
                % response time out!
                data.Response(tt) = 0;
                keepChecking=0;
            end    
                       
        end

        TimeUpdate = TimeUpdate +t.ITI(tt); %Update Matlab on what time it is.
        
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
    TimeUpdate = TimeUpdate + t.EndFixation;
    
    t.EndTime = GetSecs; %Get endtime of the experiment in seconds
    t.TotalExpTime = (t.EndTime-t.StartTime); %Gets the duration of the total run.
    t.TotalExpTimeMins = t.TotalExpTime/60; %TOTAL exp time in mins including begin and end fixation.
    t.GlobalTimer = GlobalTimer; %Spits out the exp time in secs excluding begin and end fixation.

    %% get accuracy
    
    % don't worry about time out trials or the first trial
    inds2use = (1:p.nTrials)'~=1 & data.Response~=0 & ~isnan(data.Response);
    
    data.Accuracy = mean(data.Response(inds2use)==p.correct_resp(inds2use));
    data.hit_rate = sum(data.Response(inds2use)==1 & p.correct_resp(inds2use)==1)/sum(p.correct_resp(inds2use)==1);
    data.false_alarm_rate = sum(data.Response(inds2use)==1 & p.correct_resp(inds2use)==2)/sum(p.correct_resp(inds2use)==2);
    data.timeout_rate = sum(data.Response==0)/p.nTrials;
    
    fprintf('\nCompleted block %d!\n',p.runNumGlobal);
    fprintf('Percent accuracy: %.2f\n', data.Accuracy*100);
    fprintf('Percent of repeats detected: %.2f\n', data.hit_rate*100);
    fprintf('Percent false alarms: %.2f\n', data.false_alarm_rate*100);
    fprintf('Percent time out trials: %.2f\n',data.timeout_rate*100);
    
    InstrText = sprintf('Block finished!\n\nHit rate: %.2f percent\n\nFalse alarm rate: %.2f percent',...
        data.hit_rate*100, data.false_alarm_rate*100);

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
%     if saveremote
%         save(p.fnsave_remote,'TheData');
%     end
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

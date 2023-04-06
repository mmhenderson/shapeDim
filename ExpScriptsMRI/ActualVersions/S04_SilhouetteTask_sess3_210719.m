   % SILHOUETTE LOCALIZER - SHAPE SPACE
% Subject views flashing checkerboards within an aperture that matches the
% area of space covered by any of the shape images.
% Spatial attention is directed to the entire aperture, subject must press
% a button when they see a dimming event anywhere on the checkerboard.

% WRITTEN FOR SCANNER LAPTOP

% parameters
    % Image size 24 degrees (fills almost entire screen)
    % Each trial is 7 seconds long, ITI jittered 2-8 seconds
    % 20 trials
    % total length:
        % 20*7 + 20*5 + 13+9 = 262 seconds
        % 328 TRs
        % 4:22   min:sec
    
try
%     
    echo off
    clear 
    close all hidden
       
    p.scannerlaptop = 1;
%     saveremote = 0;
    
    expdir = pwd;
    filesepinds = find(expdir==filesep);
    root = expdir(1:filesepinds(end)-1);

    % set up paths for saving
    if p.scannerlaptop || IsWindows  
        % if scanner laptop or debugging, save to a subfolder here.
        datadir_local = fullfile(root, 'Data');
%         datadir_remote = datadir_local;
    else
        % linux, behavior rooms - save to a different folder.
        datadir_local = '/home/pclexp/Documents/Maggie/shapeDim/Data/';
%          datadir_remote = '/mnt/pclexp/Maggie/shapeDim/Data/';
    end

    %% Collect information about the subject, the date, etc.
       
    % make a dialog box to enter things
    prompt = {'Debug mode?','Subject Initials','Subject Number','Contrast decrement (0-1)',...
         'Run Number','Image Set','Low-contrast background?','Low-contrast checkerboard?'};
    dlgtitle = 'Enter Run Parameters';
    dims = [1 35];
    definput = {'0','CA ','4','0.5','1','3','1','1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    p.debug = str2double(answer{1});
    p.Subject = answer{2};
    p.SubNum = sprintf('%02d',str2double(answer{3}));
    DimBy = str2double(answer{4});
    p.runNumGlobal = str2double(answer{5});
    p.ImageSet = str2double(answer{6});
    p.LowContrastBackground = boolean(str2double(answer{7}));
    p.LowContrastChecker = boolean(str2double(answer{8}));
    
    % check values of these entries
    if DimBy<0 || DimBy>1
        error('Contrast decrement must be between 0 and 1'); end
    if ~ismember(p.ImageSet,[3]); error('image set must be 3'); end

    rng('default')
    t.MySeed = sum(100*clock); 
    rng(t.MySeed);

    p.contrastDimmed = 1-DimBy; 
    
    p.expName = 'SilhouetteLocalizer';
%     p.imagedir = fullfile(root,'Stimuli',sprintf('AmpGrid%d_grey',p.ImageSet));
    p.imagedir=fullfile(root, sprintf('Stimuli/AmpGrid%d_adj_full_grey/',p.ImageSet));
    
    %% initialize my data file

    % save a copy of the currently running script here (in text format) so
    % we can look at it later if issues arise
    p.MyScript = fileread([mfilename('fullpath'),'.m']);
    
    % make sure my data dir exists
    if ~exist(datadir_local,'dir')
        mkdir(datadir_local)
    end

    t.TheDate = datestr(now,'yymmdd'); %Collect todays date (in t.)
    t.TimeStamp = datestr(now,'HHMM'); %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
    
    p.fnsave_local = fullfile(datadir_local, ['S', p.SubNum, '_' p.expName '_' t.TheDate '.mat']);

    if exist(p.fnsave_local,'file')
        load(p.fnsave_local);
        if p.runNumGlobal~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.runNumGlobal~=1            
        error('No data exists yet for this subject, check your run number')
    end
    
   
    %% set up my screen 
    
    InitializeMatlabOpenGL;  
    PsychImaging('PrepareConfiguration');
    Screen('Preference', 'SkipSyncTests', 0);
    AssertOpenGL; % bail if current version of PTB does not use
    PsychJavaTrouble;
        
    % this is about middle gray
    % contrast (black/white luminances) are defined relative to this
    % background luminance
%     p.backColor = 127;  
    if p.LowContrastBackground
        p.backColor = 50;
        p.fixColor = [0.05, 0.05, 0.05  ] * 255;
    else
        p.backColor = 127; 
        p.fixColor = [0.80,0.80,0.80]*255;
    end

    p.windowed = 0;
    s=max(Screen('Screens'));
    p.black = BlackIndex(s);
    p.white = WhiteIndex(s);
    if p.LowContrastChecker
%         p.white = 125;
        p.white = 100;    
        p.black = 0;
    end
    % Open a screen
    Screen('Preference','VBLTimestampingMode',-1);  % for the moment, must disable high-precision timer on Win apps
    multiSample=0;
    if p.windowed
        Screen('Preference', 'SkipSyncTest  s', 1);
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

    if p.scannerlaptop
        p.refreshRate = 60;
        % INNER BORE screen
        p.vDistCM = 47;
        p.screenHeightCM = 16;
%         % STAND UP SCREEN
%         p.screenHeightCM = 90; % in cm
%         p.vDistCM = 370; % in  cm
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
        % set the priority up   way high to discourage interruptions
        Priority(MaxPriority(w));
    end

    %% size/stimulus information
    %get center of screen
    p.centerPix = [(p.sRect(3) - p.sRect(1))/2, (p.sRect(4) - p.sRect(2))/2];
    p.fixSizeDeg = .2;
    p.apertureSizeDeg = 0.8 ;
    p.checkPeriodDeg = 2;
%      Deg = 14;
    p.stimHeightDeg = 24;   
    % note that this is the size of the square image that is drawn, including its grey background. 
    % The shape itself is about 2/3 of that size. So, the background is
    % drawn in a frame slightly bigger than the screen, but all the shape pixels are
    % within the bounds of the screen.
    
    % convert from degrees to pixel units
    p = deg2pix(p);
    p.fixSizePix = ceil(p.fixSizePix);
    p.apertureSizePix = ceil(p.apertureSizePix);
    
    %% timing information
    p.nTrials = 20;
    t.FlickerRateHz = 5;
    t.StimTimeTotal= 7;     
    
    t.FrameDur = 1/t.FlickerRateHz/2;   % two frames per cycle (on/off)
%     % how long do they have to respond from image onset?
%     t.RespTimeRange = 2;
%     
    t.ITIrange = [2,8];
    itis = linspace(t.ITIrange(1),t.ITIrange(2),p.nTrials);
    t.ITI = itis(randperm(length(itis)))';
   
    if p.scannerlaptop
        if p.debug
            t.StartFixation = 1;
            t.EndFixation = 1;
        else
            t.StartFixation = 13;
            t.EndFixation = 9 ;
        end
    else
        t.StartFixation = 2;
        t.EndFixation = 0;
    end
    
    % decide how often to create my dimming events
%     p.eventFreq = 0.02;
    p.eventFreq = 0.10; 
    t.nCyclesPerTrial = t.FlickerRateHz*t.StimTimeTotal;
    
    % this is the number of frames that have to elapse before another event
    % can occur (this number/2 is the number of cycles that have to elapse)
    t.MinFramesBetweenEvts = 8; 
    t.MinCyclesBetweenEvts = t.MinFramesBetweenEvts/2;
    % responses also have to happen within this relatively short window
    % (i.e can't happen after another event might have occured)
    % make this a little bit shorter to ensure we get all the misses
    t.MaxRespTime = t.MinFramesBetweenEvts*t.FrameDur - 0.01;  
    
    % onset of the first frame of each cycle (just to check my timing)
    t.stim_flips = nan(p.nTrials,t.nCyclesPerTrial);  
    
    %% Create my stimuli

    FlushEvents('keyDown');
    Screen(w,'TextFont','Helvetica');
    
    InstrText = sprintf('Loading/making stimuli. This is slow, please wait...\n\n\n\n');
    DrawFormattedText(w, InstrText, 'center', 'center', p.fixColor);
  
    Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
    Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    
    % load the silhouette mask 
    imfn = fullfile(p.imagedir, 'Silhouette_any.png');
    if exist(imfn,'file')
        im=imread(imfn);
    else
        error('image file %s not found!',imfn)
    end
    
    % set up a frame to plot the image in
    p.stimWidthPix = p.stimHeightPix*size(im,2)/size(im,1);
    p.framePos=[p.centerPix(1)-p.stimWidthPix/2,p.centerPix(2)-p.stimHeightPix/2,p.centerPix(1)+p.stimWidthPix/2,p.centerPix(2)+p.stimHeightPix/2];
 
    % initialize my patches, starting with a meshgrid
    [xgrid, ygrid] = meshgrid(1:ceil(p.stimWidthPix),1:ceil(p.stimHeightPix));
    % resize the mask image to fit in here
    mask = imresize(im, size(xgrid));
    
    all_frames = [];
    
    % how many cycles per pixel? Spatial frequency
    freq_cycles_per_pix = 1./p.checkPeriodPix;
    % how many cycles in 2pi pixels? This is the frequency of the sine
    % wave we want to draw.
    freq_cycles_per_2pi = 2*pi*freq_cycles_per_pix;
    
    % looping over trials and making texture array 
    for tt=1:p.nTrials
        % looping over cycles within each trial
        % keep track of which cycle the last dimming evt occured at
        last_dim = 1-t.MinCyclesBetweenEvts;
        for cc=1:t.nCyclesPerTrial
              
            % decide whether this is a frame I should dim
            do_dim = 0;
            if cc - last_dim >= t.MinCyclesBetweenEvts
                do_dim = rand<p.eventFreq;
                last_dim = cc; 
            end

            % choose a random phase, spanning from 0-pi
            rand_phase_rads = datasample(1:180,1)*pi/180;
            % draw the checkerboard
            stim = sign(sin(xgrid*freq_cycles_per_2pi - rand_phase_rads) .* sin(ygrid*freq_cycles_per_2pi - rand_phase_rads));
            
            % dim it if needed
            if do_dim
                stim(stim==-1) = (p.black-p.backColor)*p.contrastDimmed+p.backColor;
                stim(stim==1) = p.backColor-(p.backColor-p.white)*p.contrastDimmed;
            else
                stim(stim==-1) = p.black;
                stim(stim==1) = p.white;
            end
            % mask the checkerboard with my shape
            stim(mask==mask(1,1)) = p.backColor;
            
            % put it in a texture I can access later
            all_frames(tt,cc).text = Screen('MakeTexture',w,stim);
            all_frames(tt,cc).has_dim = do_dim;
        end
        
    end
    
    p.has_dim = reshape([all_frames.has_dim],size(all_frames));
    for tt=1:size(p.has_dim,1)
        evts = find(p.has_dim(tt,:));
        assert(all(abs(diff(evts))>=t.MinCyclesBetweenEvts));
    end
    %% keys
    KbName('UnifyKeyNames')

    %use number pad - change this for scanner 
    if p.scannerlaptop
%         p.keys=[KbName('b'),KbName('y'),KbName('g'),KbName('r')];
        p.keys=[KbName('b')];
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

    %% Allocate arrays to store trial info
    
    % these are each nTrials x nFramesTotal (2 frames per cycle, plus a
    % frame for anything that happens during the ITI)
    data.Response = nan(p.nTrials, t.nCyclesPerTrial*2+1);    % what button did they press on this trial?
    data.ResponseTimes = nan(p.nTrials, t.nCyclesPerTrial*2+1);    % when did they press the button?
    data.RespGood = nan(p.nTrials, t.nCyclesPerTrial*2+1);    % was this a good time to press or a false alarm?
    data.RespTimeFromOnset = nan(p.nTrials, t.nCyclesPerTrial*2+1);    % when was this press relative to last event (if a hit?) 
    
    % count things, just to double check responses
    data.MissCount = 0;
    data.HitCount = 0;
    data.EvtCount = 0;
    data.FrameCount = 0;
    data.FACount = 0;

    %% START EXPERIMENT
    % Draw an instruction screen, wait for space press
    FlushEvents('keyDown');
    Screen(w,'TextFont','Helvetica');
    
    InstrText = sprintf('Fixate.\nPress 1 when checkers get dimmer.\n\n\n\n');
    DrawFormattedText(w, InstrText, 'center', 'center', p.fixColor);
  
    Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
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
    Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
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
    
    % at first, no events have happened so they shouldn't be pressing
    % a button.
    keepChecking = 0;
    last_evt = nan;

    % make sure that we don't double-count any button presses
    currently_pressing = 0;

    for tt=1:p.nTrials

        
        %% Show target image   

        fc= 0;  % fc = frame count, counting all my frames (including both on/off)
        for cc=1:t.nCyclesPerTrial
  
            fc = fc+1;
            data.FrameCount = data.FrameCount+1;
            %% DRAW THE CHECKERBOARD
            Screen('DrawTexture', w, all_frames(tt,cc).text,[],p.framePos);
            Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);

            t.stim_flips(tt,cc) = GetSecs;

            Screen('Flip', w); 
            GlobalTimer = GlobalTimer + t.FrameDur;
            
            % if there is a dim here, record its time and start looking for
            % responses
            if all_frames(tt,cc).has_dim==1
                assert(keepChecking==0 && isnan(last_evt))
                keepChecking = 1;
                last_evt = t.stim_flips(tt,cc);
                data.EvtCount = data.EvtCount + 1;
            end

            TimePassed = (GetSecs-TimeUpdate);
            while TimePassed < t.FrameDur
                TimePassed = (GetSecs-TimeUpdate);
                %check for escape responses 
                [resp, timeStamp] = checkForResp(p.keys,p.escape);
                if resp==-1; escaperesponse(OriginalCLUT); end
                % only look for responses if they haven't yet responded on
                % this frame, and if they're not already holding the button
                % down
                if ~currently_pressing && isnan(data.Response(tt,fc)) && resp && find(p.keys==resp) 
                    
                    % they have responded, record it
                    data.Response(tt,fc) = find(p.keys==resp);
                    data.ResponseTimes(tt,fc) = timeStamp;
                    currently_pressing = 1; % don't look for another press until they lift the button
                    
                    % now determine if it's a false alarm or a hit
                    if keepChecking
                        % hit
                        data.HitCount = data.HitCount + 1;
                        data.RespTimeFromOnset(tt,fc) = timeStamp - last_evt;
                        data.RespGood(tt,fc) = 1;
                        % if they have responded, they should not respond
                        % again.
                        keepChecking = 0;
                        last_evt = nan;
                        
                    else
                        % false alarm
                        data.RespTimeFromOnset(tt,fc) = nan;
                        data.RespGood(tt,fc) = 0;
                        data.FACount = data.FACount + 1;
                    end
                   
                elseif currently_pressing && resp==0
                    % they were pressing the button before, but now they've
                    % lifted it, so start looking for new presses again.
                    currently_pressing = 0;
                end
                
                % finally, check if it's been too long since the last event
                % and they have missed a response
                if keepChecking && GetSecs-last_evt>=t.MaxRespTime
                    keepChecking = 0;
                    last_evt = nan;
                    data.MissCount = data.MissCount +1;
                end
            end
            TimeUpdate = TimeUpdate + t.FrameDur;

            %% DRAW THE BLANK IN-BETWEEN FRAME
            fc = fc+1;
            data.FrameCount = data.FrameCount + 1;
%             Screen('DrawTexture', w, all_frames(all_frames(tt,cc).text,[],p.framePos);
            Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);


            Screen('Flip', w); 
            GlobalTimer = GlobalTimer + t.FrameDur;

            TimePassed = (GetSecs-TimeUpdate);
            while TimePassed < t.FrameDur
                TimePassed = (GetSecs-TimeUpdate);
                %check for escape responses 
                [resp, timeStamp] = checkForResp(p.keys,p.escape);
                if resp==-1; escaperesponse(OriginalCLUT); end
                if ~currently_pressing && isnan(data.Response(tt,fc)) && resp && find(p.keys==resp) 
                    
                    % they have responded, record it
                    data.Response(tt,fc) = find(p.keys==resp);
                    data.ResponseTimes(tt,fc) = timeStamp;
                    currently_pressing = 1; % don't look for another press until they lift the button
                    
                    % now determine if it's a false alarm or a hit
                    if keepChecking
                        % hit
                        data.HitCount = data.HitCount + 1;
                        data.RespTimeFromOnset(tt,fc) = timeStamp - last_evt;
                        data.RespGood(tt,fc) = 1;
                        % if they have responded, they should not respond
                        % again.
                        keepChecking = 0;
                        last_evt = nan;
                        
                    else
                        % false alarm
                        data.RespTimeFromOnset(tt,fc) = nan;
                        data.RespGood(tt,fc) = 0;
                        data.FACount = data.FACount + 1;
                    end
                   
                elseif currently_pressing && resp==0
                    % they were pressing the button before, but now they've
                    % lifted it, so start looking for new presses again.
                    currently_pressing = 0;
                end
                
                % finally, check if it's been too long since the last event
                % and they have missed a response
                if keepChecking && GetSecs-last_evt>=t.MaxRespTime
                    keepChecking = 0;
                    last_evt = nan;
                    data.MissCount = data.MissCount +1;
                end
            end
            
            TimeUpdate = TimeUpdate + t.FrameDur;

        end
        
%         TimeUpdate = TimeUpdate + t.StimTimeTotal; %Update Matlab on what time it is.

        %% DO THE ITI - CHECK RESPONSES INTO THIS PERIOD IF NECESSARY
        Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
        Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
        Screen('DrawingFinished', w);                     
        Screen('Flip', w);
       
        %TIMING!:
        fc = fc + 1;    % consider the ITI to be the last "Frame" of the trial
        data.FrameCount = data.FrameCount + 1;
        GlobalTimer = GlobalTimer + t.ITI(tt);
        TimePassed = (GetSecs-TimeUpdate);
        while TimePassed < t.ITI(tt)
            TimePassed = (GetSecs-TimeUpdate);
            %check for escape responses
            [resp, timeStamp] = checkForResp(p.keys,p.escape);
            if resp==-1; escaperesponse(OriginalCLUT); end
            if ~currently_pressing && isnan(data.Response(tt,fc)) && resp && find(p.keys==resp)
                
                % they have responded, record it
                data.Response(tt,fc) = find(p.keys==resp);
                data.ResponseTimes(tt,fc) = timeStamp;
                currently_pressing = 1; % don't look for another press until they lift the button
                
                % now determine if it's a false alarm or a h it
                if keepChecking
                    % hit
                    data.HitCount = data.HitCount + 1;
                    data.RespTimeFromOnset(tt,fc) = timeStamp - last_evt;
                    data.RespGood(tt,fc) = 1;
                    % if they have responded, they should not respond
                    % again.
                    keepChecking = 0;
                    last_evt = nan;
                    
                else
                    % false alarm
                    data.RespTimeFromOnset(tt,fc) = nan;
                    data.RespGood(tt,fc) = 0;
                    data.FACount = data.FACount+1;
                end
                
            elseif currently_pressing && resp==0
                % they were pressing the button before, but now they've
                % lifted it, so start looking for new presses again.
                currently_pressing = 0;
            end
            
            % finally, check if it's been too long since the last event
            % and they have missed a response
            if keepChecking && GetSecs-last_evt>=t.MaxRespTime
                keepChecking = 0;
                last_evt = nan;
                data.MissCount = data.MissCount +1;
            end
        end

        TimeUpdate = TimeUpdate +t.ITI(tt); %Update Matlab on what time it is.
        
    end 
    
    %% finish experiment 
       
    % final fixation:
    Screen('DrawDots', w, [0,0], p.apertureSizePix, p.backColor, p.centerPix, 1);
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
    
    % check total events 
    total_evts = sum(p.has_dim(:));
    assert(total_evts==data.EvtCount);
    assert(data.EvtCount==data.HitCount+data.MissCount);
     
    % check false alarms
    false_alarms = sum(data.RespGood(:)==0);
    assert(false_alarms==data.FACount);
    
    % check hits
    hits_per_trial = sum(data.RespGood==1,2);
    assert(sum(hits_per_trial)==data.HitCount);
    data.hit_rate = sum(hits_per_trial)/total_evts;
    
    % check misses
    misses_per_trial = sum(p.has_dim,2) - sum(data.RespGood==1,2);
    assert(sum(misses_per_trial)==data.MissCount);
    data.miss_rate = sum(misses_per_trial)/total_evts;
    
    % looking at miss trials and seeing how late they were
%     tt = 6;
%     evt_times = t.stim_flips(tt,p.has_dim(tt,:)==1);
%     resp_times = data.ResponseTimes(tt,data.Response(tt,:)==1);
%     resp_offset = resp_times-evt_times;
    
    fprintf('\nCompleted block %d!\n',p.runNumGlobal);
    fprintf('Hit rate: %.2f percent',data.hit_rate*100);    
    fprintf('Number false alarms: %d\n', data.FACount);
    
    InstrText = sprintf('Block finished!\n\nHit rate: %.2f percent\n\nFalse alarms: %d',...
        data.hit_rate*100, data.FACount);

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
    %
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
   
    sca
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


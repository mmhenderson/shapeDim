 % CATEGORY JUDGMENT TASK - SHAPE SPACE     
% Subject sees a sequence of images - must report the category of each
% image, according to a two-way division. 
% There are several divisions that can be used.
% We'll present two prototype images to remind the subject which images are
% which.   
 
% WRITTEN FOR BEHAVIOR ROOMS AND SCANNER LAPTOP
% change variable p.scannerlaptop to indicate which situation.
        
% parameters for MRI task
    % Image size 24 degrees (fills almost entire screen)
    % Image onscreen for 1 second    
    % Response period for 1 second (blank)
    % ITI [1-5] seconds jittered
    
    % 48 total images presented
    
    % total time = 13s blank start + 48*(1+1+3) + 8s blank end 
    %   = 261 seconds or 4:21
    %   with TR=0.8 sec this is 327 TRs.
    
try
    
    echo off
    clear 
    close all hidden
       
    p.scannerlaptop = 1;

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
        'Part (1-3)', 'Run Number (1-3)', 'Training Run? (0 or 1)','Image Set','Difficulty (7-13)'};
    dlgtitle = 'Enter Run Parameters';
    dims = [1 35];
    definput = {'0','CP','1','3','1','1','0','3','9'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    p.debug = str2double(answer{1});
    p.SubInit = answer{2};  
    p.SubNum = str2double(answer{3});
    p.SubNumStr = sprintf('%02d',str2double(answer{3}));
    p.Session = str2double(answer{4});
    p.Part = str2double(answer{5});
    p.RunNumThisPart = str2double(answer{6});
    p.Training = str2double(answer{7});
    p.ImageSet = str2double(answer{8});
    p.RunDifficulty = str2double(answer{9});
    
    % check values of these entries
    if ~ismember(p.Session,[1:3]);  error('session must be 1, 2, or 3');  end
    if ~ismember(p.Part,1:6); error('part must be 1-6'); end
    if ~ismember(p.Training,[0,1]); error('training must be 0 or  b 1'); end
    if ~ismember(p.ImageSet,[3]); error('image set must be 3'); end
    if ~ismember(p.RunDifficulty,7:13); error('difficulty level must be 7-13'); end
   
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
         
    if p.Training
        p.expName = sprintf('MainTaskMRI_scannerversion_sess%d_part%d_TRAINING',p.Session,p.Part);
    else
        p.expName = sprintf('MainTaskMRI_scannerversion_sess%d_part%d',p.Session,p.Part);
    end
    
    p.fnsave_local = fullfile(datadir_local, ['S', p.SubNumStr, '_' p.expName '_' t.TheDate '.mat']);
    p.fnsave_remote = fullfile(datadir_remote, ['S', p.SubNumStr, '_' p.expName '_' t.TheDate '.mat']);
    
    if exist(p.fnsave_local,'file')
        load(p.fnsave_local);
        if p.RunNumThisPart~=length(TheData)+1
             error('Check your run number, %d runs have been completed',length(TheData));
        end
    elseif p.RunNumThisPart~=1            
        error('No data exists yet for this subject, check your run number')
    end
    
    p.rndseed = round(sum(100*clock));
    rng(p.rndseed);
    
    if p.debug || p.Training
        % very small ## trials here
        p.nTrials = 5;
        p.nTrials_variable = 0;
        p.nTrials_main = 5;
    else   
        p.nTrials = 48;
        p.nTrials_variable = 16;    % variable are the trials whose difficulty is set by experimenter
        p.nTrials_main = 32;    % main are trials from the fixed 16-pt grid.
    end

    %% figure out what to do during this block.
    
    % list the orders in which things will happen across entire experiment
    % [nParts x nSess]
    
    % which boundary is active? this defines which categorization dim they're using
    bound_list = [1,2,3; 2,3,1; 3,1,2; 1,2,3; 2,3,1; 3,1,2];
    % which response mapping? e.g. which finger for which category?
    map_list = [1,2,1; 1,2,2; 1,2,1; 2,1,2; 2,1,1; 2,1,2];

    if mod(p.SubNum,3)==2
        % second session goes first
        bound_list = bound_list(:,[2,3,1]);
        map_list = map_list(:,[2,3,1]);
    elseif mod(p.SubNum,3)==0 
        % third session goes first
        bound_list = bound_list(:,[3,1,2]);
        map_list = map_list(:,[3,1,2]);
    end
    % else keep the orig order, first session first 
        
    % check this order to make sure we are including all 6 tasks in the current session
    unrows = unique([bound_list(:,p.Session),map_list(:,p.Session)],'rows');
    assert(size(unrows,1)==6);
    position_list = map_list;
    
    % choose the right parameters for this part
    p.which_bound = bound_list(p.Part, p.Session);
    p.which_mapping = map_list(p.Part, p.Session);
    p.which_pos = position_list(p.Part, p.Session);
    
    %% set up the main 16 pts grid 
   
    % fixed number of images in the grid
    p.nIms_main = 16;
    % Define a 4x4 grid w even spacing
    start = 0;    % min value along each axis
    stop = 5; % max value along each axis  
    nsteps_main = 4;    % how many steps along each axis?
    start_grid = 0.1;
    stop_grid = 4.9;
    main_pts = round(linspace(start_grid,stop_grid, nsteps_main),1);
    
    [gridx,gridy] = meshgrid(main_pts,main_pts);
    main_grid_points = [gridx(:),gridy(:)];
    p.main_grid_points = main_grid_points;
    
    %% set up the full image grid 
   
    % first I'm defining all possible images that we can use in this task. 
    step = 0.1;
    center = (stop-start)./2+start;
      
    all_pts = round(start:step:stop,1) ;  
    [gridx,gridy] = meshgrid(all_pts,all_pts);
    all_grid_points = [gridx(:),gridy(:)];
    
    % now taking out images at exactly the prototype locations so that we
    % never use these during task
    proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
    proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
    proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
    assert(numel(proto_inds)==4)
    all_grid_points(proto_inds,:) = [];
  
    % Also taking out any images along quadrant boundaries because these
    % are ambiguous 
    bound_inds = find(any(all_grid_points==center,2));
    all_grid_points(bound_inds,:) = [];
  
    % now define which quadrant each image lies in. note that this
    % "quadrant" property is fixed no matter what the task is, the
    % "category" is a separate property that will change depending on task.
    % that gets defined later.
    all_quadrant = zeros(size(all_grid_points,1),1);
    all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)>center) = 1;
    all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)>center) = 2;
    all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)<center) = 3;
    all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)<center) = 4;
    
    % Next, for each point in the full grid, define the difficulty level 
    % based on distance from the boundary. This depends on the active task
    % because it's relative to the active boundary.
    if p.which_bound<3
        % just considering distance from the binary bound
        dist_from_bound = round(abs(all_grid_points(:,p.which_bound)-center),1);
    else
        % considering distance from any boundary
        dist_from_bound = round(min(abs((all_grid_points-repmat(center, size(all_grid_points,1),2))),[],2),1);
    end
    
    % bin these values into 13 "difficulty levels"
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
    assert(all(nperbin>=p.nTrials_variable));
    
    %% Now choose the images that we want to show on this block. 

    % first, the 16 main grid points
    main_grid_inds = zeros(p.nIms_main,1);
    for pp=1:size(p.main_grid_points,1)
       main_grid_inds(pp) = find(ismember(all_grid_points, p.main_grid_points(pp,:),'rows'));    
    end
    main_grid_dist = dist_groups_binned(main_grid_inds);
    
    % next the variable trials
    % based on the difficulty level inputted above, decide what shape space
    % positions to sample.
    difficulty_to_sample = [p.RunDifficulty];
    p.difficulty_to_sample = difficulty_to_sample;
    
    assert(~mod(p.nTrials_variable, length(difficulty_to_sample)));
    assert(~mod(p.nTrials_variable, 4));
    nEach = p.nTrials_variable/length(difficulty_to_sample)/4;
    var_inds = [];
    for dd=1:length(difficulty_to_sample)
        % looping over quadrants to keep the responses balanced
        for qq=1:4           
            % from each quadrant and each difficulty level that we want to
            % sample, grab a random subset of image indices.
            thisbin = find(dist_groups_binned==difficulty_to_sample(dd) & all_quadrant==qq);    
            var_inds = [var_inds; datasample(thisbin, nEach, 'replace',false)];
        end
    end
  
    % check this math
    points_check = all_grid_points(var_inds,:); 
    if p.which_bound<3
        dist_from_bound_check = round(abs(points_check(:,p.which_bound)-center),1);
    else
        dist_from_bound_check = round(min(abs((points_check-repmat(center, size(points_check,1),2))),[],2),1);
    end
    assert(all(dist_from_bound_check==repmat(bin_edges(difficulty_to_sample)',nEach*4,1))); 
        
    % now combine the main grid images and the variable diff images
    if p.Training==1 || p.debug==1
        myinds = [datasample(main_grid_inds,p.nTrials_main)];
    else
        % for full version of task, using each of main images 2x
        myinds = [main_grid_inds; main_grid_inds; var_inds];  
    end
    is_main_grid = [ones(p.nTrials_main,1);zeros(p.nTrials_variable,1)];
    assert(numel(myinds)==p.nTrials);
    assert(numel(is_main_grid)==p.nTrials);
    
    % index into the big lists to get all the shapes we want to use and
    % their properties
    p.points = all_grid_points(myinds,:);
    p.dist_from_bound = dist_from_bound(myinds);
    p.trial_difficulty = dist_groups_binned(myinds);
    p.quadrant = all_quadrant(myinds);
    p.is_main_grid = is_main_grid;
        
    %% Make a randomized image sequence    
    trial_order = (1:p.nTrials)';
    trial_order = trial_order(randperm(length(trial_order)));
    
    % sort everything the same way
    p.points = p.points(trial_order,:);
    p.dist_from_bound= p.dist_from_bound(trial_order);
    p.trial_difficulty = p.trial_difficulty(trial_order);
    p.quadrant = p.quadrant(trial_order);
    p.is_main_grid = p.is_main_grid(trial_order);
      
    %% Define the categories/response mapping  for this run
    % which categorization scheme are we using? each row defines the two
    % "quadrants" that are grouped into a common "category"
    if p.which_bound == 1
        cat_groups = [1, 4; 2, 3];
    elseif p.which_bound == 2
        cat_groups = [1, 2; 3, 4];
    else
        cat_groups = [1, 3; 2, 4];
    end        
    
    if p.which_mapping == 2
        % flip the rows, will flip which finger corresponds to which
        % category. At start of task, the prototypes for the categories
        % will always be shown on L/R side of fixation, so changing the
        % mapping will also switch which side of fixation the categories
        % are each shown on.
        cat_groups = cat_groups([2,1],:);
    end
    
    if p.which_pos == 2
        % flip the spatial positions in which the prototypes are shown
        % (within their groups). Note that this is NOT counter-balanced
        % w/r/t which_mapping; they always have the same value in this vers
        % of the code. Left over from older vers of this script, not super
        % important. 
        cat_groups = cat_groups(:,[2,1]);
    end
    
    % define positions to draw prototypes on start screen.
    p.proto_order = [cat_groups(1,:), cat_groups(2,:)];
    
    % define category for each trial. 
    % quadrant is a fixed property of each image
    % category depends on the rule in place (which finger/side)
    p.category = zeros(size(p.quadrant));
    p.category(ismember(p.quadrant, cat_groups(1,:)))= 1;
    p.category(ismember(p.quadrant, cat_groups(2,:)))= 2;
    
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
%     line([center,center],get(gca,'YLim'),'Color','k')
%     line(get(gca,'XLim'), [center,center], 'Color','k');
%     
%     set(gcf,'Color','w');
%     title(sprintf('Shape locations used for this block: Boundary %d, Map %d, Difficulty %d\n',p.which_bound, p.which_mapping, p.RunDifficulty));
    
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
    p.fixColor = [0.8,0.8,0.8]*255;
%     p.fixColor = [0,0,0];
    p.stimHeightDeg = 24;   
    % note that this is the size of the square image that is drawn, including its grey background. 
    % The shape itself is about 2/3 of that size. So, the background is
    % drawn in a frame slightly bigger than the screen, but all the shape pixels are
    % within the bounds of the screen.

    if  p.scannerlaptop
        % have to shrink the prototypes a tiny bit to accomodate the screen
        % size here 
        p.protoSpacingDeg=0.01;
        p.protoStimHeightDeg=6 ;
%         p.protoSpacingDeg = 0.01; % how far apart are the four prototypes drawn?
%         p.protoStimHeightDeg = 8.5;  % how big is each prototype image?    
        p.feedbackTextHeightDeg = 5 ;
    else
        p.protoSpacingDeg = 0.01; % how far apart are the four prototypes drawn?
        p.protoStimHeightDeg = 8.5;  % how big is each prototype image?    
        p.feedbackTextHeightDeg = 9 ;
    end
    p.outlineProp = 0.9; % for the colored square outlines, how big are they relative to the image square size?
    p.rectColor = [0.1,0.8,0.1]*255;
    p.rectWidthDeg = 0.1;
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
%     protocenters_pix = repmat(p.centerPix,4, 1) + p.protoDistPix*[-3,0; -1.2,0; 1.2,0; 3 ,0];    
   
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
    % now create the two boxes, each spanning two images.
    p.catOutlinePos = [outlines(1,1),outlines(1,2),outlines(2,3),outlines(1,4);...
                                    outlines(3,1),outlines(3,2),outlines(4,3),outlines(3,4)];

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
    t.ShowFeedbackTime = 1;
%     
    if p.Training
        t.MaxRespTime = 10;  
    else
        % how long do they have to respond from image onset?
        t.RespTimeRange = 2;
        t.RespTimeBlank = t.RespTimeRange - t.StimTimeTotal;
    end
    
    if p.scannerlaptop && p.Training==0
        % actual timing for scanner expt - long delay at start,
        % longer/jittered ITIs
        t.StartFixation = 13;
        t.EndFixation = 8;
        t.ITIrange = [1, 5];
    else
        % for behav room, no need to have these longer delays so whole expt
        % will be shorter. Same thing for mini training blocks at scanner,
        % bc scanner not actually on for these.
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
 
    % draw the prototype images
    for bb=1:4
        Screen('DrawTexture', w, proto_ims(p.proto_order(bb)).imtext,[],p.protoFramePos(bb,:));        
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
       
        % next part depends on whether training or real block...
        if p.Training
            
            %% draw the prototype images on screen
            % to remind subject of response mapping.
            for bb=1:4
                Screen('DrawTexture', w, proto_ims(p.proto_order(bb)).imtext,[],p.protoFramePos(bb,:));        
%                 DrawFormattedText(w, sprintf('%d',bb), p.boundTextPos{oo}(bb,1)+p.textAdjPix(1),p.boundTextPos{oo}(bb,2)+p.textAdjPix(2),p.white);
            end

            Screen('DrawDots', w, [0,0], p.fixSizePix, p.fixColor, p.centerPix, 0); 
            Screen('DrawingFinished', w);

            t.stim_flips(tt,2) = GetSecs;

            Screen('Flip', w);
            
            % will give a maximum time window to respond, if they respond
            % before then they move onto next trial 
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
                Screen('DrawTexture', w, proto_ims(p.proto_order(bb)).imtext,[],p.protoFramePos(bb,:));        
%                 DrawFormattedText(w, sprintf('%d',bb), p.boundTextPos{oo}(bb,1)+p.textAdjPix(1),p.boundTextPos{oo}(bb,2)+p.textAdjPix(2),p.white);
            end
                       
            if data.Response(tt)~=0
                 % outline their answer
                Screen('FrameRect', w, p.fixColor, p.catOutlinePos(data.Response(tt),:),p.rectWidthPix);
             
                if data.Response(tt)==p.category(tt)
                    showtext = 'Correct!\n\n';
                else
                    showtext = 'Incorrect.\n\nCorrect response is shown in green.';
                end
            else 
               showtext = 'Time out!\n\nCorrect response is shown in green.'; 
            end
            
            % outline the correct answer
            Screen('FrameRect', w, p.rectColor, p.catOutlinePos(p.category(tt),:),p.rectWidthPix);
                          
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
                    
        else
            %% FULL TASK - no reminders of prototypes
            % just blank screen, but keep checking responses

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

        end
        
        % now back to same for training/main...
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
    
    data.MainGridAccuracy = mean(data.Response(trialsdone & p.is_main_grid==1)==...
        p.category(trialsdone & p.is_main_grid==1));
    data.VariableAccuracy = mean(data.Response(trialsdone & p.is_main_grid==0)==...
        p.category(trialsdone & p.is_main_grid==0));
    
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
    t.TotalExpTime = (t.EndTime-t.StartTime); %Gets the duration of the total run
    t.TotalExpTimeMins = t.TotalExpTime/60; %TOTAL exp time in mins including begin and end fixation.
    t.GlobalTimer = GlobalTimer; %Spits out the exp time in secs excluding begin and end fixation.

    %% get accuracy
    
    fprintf('\nCompleted block %d!\n',p.RunNumThisPart);
    fprintf('Accuracy is %.2f percent\n',data.Accuracy*100);
    fprintf('Easier trials: %.2f percent\n',data.MainGridAccuracy*100);
    fprintf('Harder trials (diff=%d): %.2f percent\n',p.RunDifficulty,data.VariableAccuracy*100);
    fprintf('Number of time out trials: %d/%d\n',sum(data.Response==0),p.nTrials);
    
    InstrText = ['Block finished!' '\n\n'...
                sprintf('Accuracy is %.2f percent',data.Accuracy*100)];
                

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
    TheData(p.RunNumThisPart).t = t;
    TheData(p.RunNumThisPart).p = p;
    TheData(p.RunNumThisPart).data = data;
 
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


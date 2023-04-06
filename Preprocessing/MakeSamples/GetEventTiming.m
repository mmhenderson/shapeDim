function [] =  GetEventTiming(sub2do)

    % Get the timing of all events in all runs
    % saves out a file that is used by MakeSampleFile and by RunGLM_WholeBrain

    % MMH 9/5/18

    clear
    close all

    if nargin<1
        sub2do=[1,2,3,4,5,6,7];
    end
    
    subinit_big = {'CP','BX','BR','CA','CG','CT','CU'};
    subnum_big = [1,2,3,4,5,6,7];
    % sub2do = [1,2,3,4,5,6,7];
    % sub2do = [6];

     % set paths (end all with filesep!)
    mypath = pwd;
    filesepinds = find(mypath==filesep);
    nDirsUp = 2;
    experiment_path = mypath(1:filesepinds(end-nDirsUp+1));

    out_path = fullfile(experiment_path,'Samples');
    beh_path = fullfile(experiment_path, 'DataBehavior');

    for ss = sub2do


        % set inputs
        FS_sbj = subinit_big(ss);
        subnum = subnum_big(ss);
        substr = sprintf('S%02d',subnum);

        % these are files that i made manually, listing which runs we need to
        % load for the main task. For some subjects, we're loading a run from
        % session 3 in place of one from session 1, etc. so this is addressed
        % in this file...
        main_run_sequence = csvread(fullfile(out_path,'run_sequences',sprintf('%s_main_run_sequence.csv',substr)),1,0);
        main_sess_num_load = main_run_sequence(:,4);
        main_sess_date_load = main_run_sequence(:,6);

        %% Timing information

        TRdur = .8;
        TRditched = 16; % 16 TRs were removed
        not_recorded = TRdur * TRditched;

        %% Main Task

        nRuns = 0; nTRs_main = 327 - TRditched; 
        RunLabels = []; EventLabels = []; TrialLabels = []; SessLabels = [];
        BoundLabels = []; MapLabels = []; 
        PointLabels = []; QuadrantLabels = []; CatLabels = [];
        IsMainLabels = [];  DistLabels = [];    DiffLabels = [];    RunDiffLabels = [];
        RTLabels = []; ResponseLabels = [];

        trialnumglobal = 0;

        nSess = 3;
        nParts = 6;
        nRunsPerPart = 2;

        % get data file from each session (they are not concatenated here)
        main_files = dir(fullfile(beh_path,substr,'Session*','*MainTask*.mat'));
        main_files = main_files(~contains({main_files.name},'TRAINING'));
        assert(length(main_files)==nSess*nParts)
        names = {main_files.name};
        dates = cellfun(@(x)x(end-9:end-4),names,'UniformOutput',false);
        sess_dates = unique(dates);

        for sess = 1:nSess

            for part  = 1:nParts

                count = nRuns+1;

                if count>length(main_sess_num_load)
                    fprintf('skipping run 1&2 of part %d for %s sess %d (no mri data)\n', part, substr, sess)
                    continue
                end

                fn2load = fullfile(beh_path,substr,sprintf('Session%d',main_sess_num_load(count)),...
                    sprintf('%s_MainTaskMRI_scannerversion_sess%d_part%d_%d.mat',substr,sess,part,main_sess_date_load(count)));
        %         if subnum==2 && strcmp(FS_sbj, 'BX') && sess==1 && (part==4 || part==6)
        %             % subject BX, we re-did a few runs of sess 1 as part of sess 3.
        %             fn2load = fullfile(beh_path,substr,sprintf('Session%d',3),...
        %                 sprintf('%s_MainTaskMRI_scannerversion_sess%d_part%d_%s.mat',substr,sess,part,sess_dates{3}));
        %         end
        %         if subnum==7 && strcmp(FS_sbj, 'CU') && ((sess==1 && part>=5) || (sess==2 && part==6))
        %              % subject CU, we re-did a few runs of sess 1/2 as part of sess 3.
        %             fn2load = fullfile(beh_path,substr,sprintf('Session%d',3),...
        %                 sprintf('%s_MainTaskMRI_scannerversion_sess%d_part%d_%s.mat',substr,sess,part,sess_dates{3}));
        % 
        %         end
                fprintf('loading from %s\n',fn2load)
                load(fn2load);

                if sess==1 && part==1
                    nTrialsEach = TheData(1).p.nTrials;
                end

                for run = 1:nRunsPerPart

                    nRuns = nRuns + 1;

                    if nRuns>length(main_sess_num_load)
                        fprintf('skipping run %d part %d for %s sess %d (no mri data)\n', run, part, substr, sess)
                        continue
                    end
    %                 if run>length(TheData)
    %                     fprintf('missing part %d run %d for %s sess %d\n', run, part, substr, sess)
    %                     continue
    %                 end


                    % Find stimulus onset time
                    Trial_onset = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime  -not_recorded;
                    Trial_offset = TheData(run).t.stim_flips(:,2) - TheData(run).t.StartTime - not_recorded;
                    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;

                    for n = 1:TheData(run).p.nTrials

                        trialnumglobal = trialnumglobal+1;

                        % list times for all events on this trial 
                        event_times_this_trial = [Trial_onset(n),Trial_offset(n)];

                        Event_onset = [Event_onset, event_times_this_trial];
                        % 1 is stim on, 0 is stim off
                        Event_type = [Event_type, 1, 0]; 
                        % also appending to a list, that is nEvents long and describes which trial each event is from 
                        Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,2)];  
                        Event_trial_local = [Event_trial_local, repmat(n,1,2)];

                    end

                    % Find onset times for all my TRs
                    TR_onset = 0:.8:nTRs_main*0.8;
                    TR_onset = TR_onset(1:end-1);

                    % list everything by TRs
                    triallabslocal_byTR = zeros(length(TR_onset),1);
                    eventlabs_byTR = zeros(length(TR_onset),1);
                    triallabsglobal_byTR = zeros(length(TR_onset),1);

                    for i = 1:length(TR_onset)

                        middleofTR = TR_onset(i)+(TRdur/2);
                        myind = find(middleofTR>Event_onset, 1, 'last');

                        % on this TR - what type of of event is happening, and what global
                        % trial number and run number is it a part of?
                        eventlabs_byTR(i) = Event_type(myind);
                        triallabslocal_byTR(i) = Event_trial_local(myind);
                        triallabsglobal_byTR(i) = Event_trial_global(myind);

                    end

                    % make sure that there is a 1 marking the start of every image. This is necessary
                    % because sometimes the events might not land in the right place in a
                    % TR
                    for n = 1:nTrialsEach           
                        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
                    end


                    % now round all the trial numbers back to being integers (need this for
                    % the code below where we use them as indices)
                    triallabslocal_byTR = round(triallabslocal_byTR,0);
                    triallabsglobal_byTR = round(triallabsglobal_byTR,0);

                    % how many TRs elapse before the first event? use this to keep all my
                    % arrays the same length
                    numZeros = sum(triallabslocal_byTR==0);

                    % append to these big lists, which are nTRs long by the end
                    TrialLabels = [TrialLabels;triallabsglobal_byTR];    
                    EventLabels = [EventLabels; eventlabs_byTR];

        %             fprintf('    boundary is %d\n',TheData(run).p.which_bound)
                    % defining some things here that are true across the entire
                    % run of this task
                    RunLabels = [RunLabels;repmat(nRuns,length(TR_onset),1)];
                    SessLabels = [SessLabels; sess*ones(size(eventlabs_byTR))];
                    BoundLabels = [BoundLabels;  TheData(run).p.which_bound*ones(size(eventlabs_byTR))];
                    MapLabels = [MapLabels;  TheData(run).p.which_mapping*ones(size(eventlabs_byTR))];
                    RunDiffLabels = [RunDiffLabels;  TheData(run).p.RunDifficulty*ones(size(eventlabs_byTR))];


                    % here i'm taking out info from my data structure, using my trial
                    % indices created above to get from trial space to TR space           
                    PointLabels = [PointLabels;  nan(numZeros,2); TheData(run).p.points(triallabslocal_byTR(numZeros+1:end),:)];
                    QuadrantLabels = [QuadrantLabels; nan(numZeros,1); TheData(run).p.quadrant(triallabslocal_byTR(numZeros+1:end))];
                    CatLabels = [CatLabels; nan(numZeros,1); TheData(run).p.category(triallabslocal_byTR(numZeros+1:end))];
                    IsMainLabels = [IsMainLabels; nan(numZeros,1); TheData(run).p.is_main_grid(triallabslocal_byTR(numZeros+1:end))];
                    DistLabels = [DistLabels; nan(numZeros,1); TheData(run).p.dist_from_bound(triallabslocal_byTR(numZeros+1:end))];
                    DiffLabels = [DiffLabels; nan(numZeros,1); TheData(run).p.trial_difficulty(triallabslocal_byTR(numZeros+1:end))];

                    RTLabels = [RTLabels; nan(numZeros,1); TheData(run).t.RespTimeFromOnset(triallabslocal_byTR(numZeros+1:end))];
                    ResponseLabels = [ResponseLabels; nan(numZeros,1); TheData(run).data.Response(triallabslocal_byTR(numZeros+1:end))];
                end

            end %run loop

        end

        % now i have a bunch of lists nTRsTotal long, and they will all get put
        % into my main data structure for saving. 

        main.RunLabels = RunLabels;
        main.TrialLabels = TrialLabels;
        main.EventLabels = EventLabels;
        main.SessLabels = SessLabels;

        main.PointLabels = PointLabels;
        main.CatLabels = CatLabels;
        main.QuadrantLabels = QuadrantLabels;

        main.BoundLabels = BoundLabels;
        main.MapLabels = MapLabels;
        main.IsMainLabels = IsMainLabels;
        main.DistLabels = DistLabels;
        main.DiffLabels = DiffLabels;
        main.RunDiffLabels = RunDiffLabels;

        main.RTLabels = RTLabels;
        main.ResponseLabels = ResponseLabels;

        if numel(main.EventLabels)~=nTRs_main*numel(unique(main.RunLabels))
            error('wrong number of total TRs!')
        end
    %     if max(main.TrialLabels)~=nSess*nParts*nRunsPerPart*nTrialsEach
    %         error('wrong number of total trials!')
    %     end

        %% Repeat detection (one-back) task

        nRuns = 0; nTRs_rep = 329 - TRditched; 
        RunLabels = []; EventLabels = []; TrialLabels = []; SessLabels = [];
        MapLabels = []; 
        PointLabels = []; IsRepeatLabels = [];
        IsMainLabels = [];  DistLabels = [];  RunDiffLabels = [];
        RTLabels = []; ResponseLabels = [];

        trialnumglobal = 0;

        nSess = 3;
        nRunsPerSess = 4;

        % get data file from each session (they are not concatenated here)
        rep_files = dir(fullfile(beh_path,substr,'Session*','*OneBackTask*.mat'));
        rep_files = rep_files(~contains({rep_files.name},'TRAINING'));
        assert(length(rep_files)==nSess)
        names = {rep_files.name};
        dates = cellfun(@(x)x(end-9:end-4),names,'UniformOutput',false);
        sess_dates = unique(dates);


        for sess = 1:nSess


            fn2load = fullfile(beh_path,substr,sprintf('Session%d',sess),...
                sprintf('%s_OneBackTaskMRI_sess%d_%s.mat',substr,sess,sess_dates{sess}));
            fprintf('loading from %s\n',fn2load)
            load(fn2load);

            for run = 1:nRunsPerSess

                if run>length(TheData)
                    fprintf('missing run %d for %s sess %d\n', run, substr, sess)
                    continue
                end


                nRuns = nRuns + 1;

                % Find stimulus onset time
                Trial_onset = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime  -not_recorded;
                Trial_offset = TheData(run).t.stim_flips(:,2) - TheData(run).t.StartTime - not_recorded;
                Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;

                for n = 1:TheData(run).p.nTrials

                    trialnumglobal = trialnumglobal+1;

                    % list times for all events on this trial 
                    event_times_this_trial = [Trial_onset(n),Trial_offset(n)];

                    Event_onset = [Event_onset, event_times_this_trial];
                    % 1 is stim on, 0 is stim off
                    Event_type = [Event_type, 1, 0]; 
                    % also appending to a list, that is nEvents long and describes which trial each event is from 
                    Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,2)];  
                    Event_trial_local = [Event_trial_local, repmat(n,1,2)];

                end

                % Find onset times for all my TRs
                TR_onset = 0:.8:nTRs_rep*0.8;
                TR_onset = TR_onset(1:end-1);

                % list everything by TRs
                triallabslocal_byTR = zeros(length(TR_onset),1);
                eventlabs_byTR = zeros(length(TR_onset),1);
                triallabsglobal_byTR = zeros(length(TR_onset),1);

                for i = 1:length(TR_onset)

                    middleofTR = TR_onset(i)+(TRdur/2);
                    myind = find(middleofTR>Event_onset, 1, 'last');

                    % on this TR - what type of of event is happening, and what global
                    % trial number and run number is it a part of?
                    eventlabs_byTR(i) = Event_type(myind);
                    triallabslocal_byTR(i) = Event_trial_local(myind);
                    triallabsglobal_byTR(i) = Event_trial_global(myind);

                end

                % make sure that there is a 1 marking the start of every image. This is necessary
                % because sometimes the events might not land in the right place in a
                % TR
                for n = 1:nTrialsEach           
                    eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
                end


                % now round all the trial numbers back to being integers (need this for
                % the code below where we use them as indices)
                triallabslocal_byTR = round(triallabslocal_byTR,0);
                triallabsglobal_byTR = round(triallabsglobal_byTR,0);

                % how many TRs elapse before the first event? use this to keep all my
                % arrays the same length
                numZeros = sum(triallabslocal_byTR==0);

                % append to these big lists, which are nTRs long by the end
                TrialLabels = [TrialLabels;triallabsglobal_byTR];    
                EventLabels = [EventLabels; eventlabs_byTR];

                % defining some things here that are true across the entire
                % run of this task
                RunLabels = [RunLabels;repmat(nRuns,length(TR_onset),1)];
                SessLabels = [SessLabels; sess*ones(size(eventlabs_byTR))];        
                MapLabels = [MapLabels;  TheData(run).p.RespMap*ones(size(eventlabs_byTR))];
                RunDiffLabels = [RunDiffLabels;  TheData(run).p.RunDifficulty*ones(size(eventlabs_byTR))];

                % here i'm taking out info from my data structure, using my trial
                % indices created above to get from trial space to TR space           
                PointLabels = [PointLabels;  nan(numZeros,2); TheData(run).p.points(triallabslocal_byTR(numZeros+1:end),:)];        
                IsRepeatLabels = [IsRepeatLabels; nan(numZeros,1); TheData(run).p.is_repeat(triallabslocal_byTR(numZeros+1:end))];
                IsMainLabels = [IsMainLabels; nan(numZeros,1); TheData(run).p.is_main_grid(triallabslocal_byTR(numZeros+1:end))];
                DistLabels = [DistLabels; nan(numZeros,1); TheData(run).p.dist_from_previous(triallabslocal_byTR(numZeros+1:end))];

                RTLabels = [RTLabels; nan(numZeros,1); TheData(run).t.RespTimeFromOnset(triallabslocal_byTR(numZeros+1:end))];
                ResponseLabels = [ResponseLabels; nan(numZeros,1); TheData(run).data.Response(triallabslocal_byTR(numZeros+1:end))];
            end

        end

        % now i have a bunch of lists nTRsTotal long, and they will all get put
        % into my main data structure for saving. 

        rep.RunLabels = RunLabels;
        rep.TrialLabels = TrialLabels;
        rep.EventLabels = EventLabels;
        rep.SessLabels = SessLabels;

        rep.PointLabels = PointLabels;
        rep.MapLabels = MapLabels;
        rep.IsRepeatLabels = IsRepeatLabels;
        rep.IsMainLabels = IsMainLabels;
        rep.DistLabels = DistLabels;
        rep.RunDiffLabels = RunDiffLabels;

        rep.RTLabels = RTLabels;
        rep.ResponseLabels = ResponseLabels;

        if numel(rep.EventLabels)~=nTRs_rep*numel(unique(rep.RunLabels))
            error('wrong number of total TRs!')
        end
        % if max(rep.TrialLabels)~=nSess*nRunsPerSess*nTrialsEach
        %     error('wrong number of total trials!')
        % end


        %% Save timing file
        filename = fullfile(out_path, sprintf('TimingFile_%s.mat',substr));
        fprintf('saving file to %s\n',filename);

        if ~exist(out_path, 'dir'), mkdir(out_path); end
        save(filename, 'main','rep');

    end

end

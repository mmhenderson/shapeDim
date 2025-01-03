%%%%%%%%%%%%%%% FIRST SCRIPT TO RUN FOR LOCALIZER ANALYSES %%%%%%%%%%%%%%%%
% get timing for each localizer run, and save out to text files that can be
% read by FEAT (one file per explanatory viariable per run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


% sublist = [1,2,3,4,5,6,7];
sublist=[8,9,10];
my_dir = pwd;

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 1;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));


for ss = 1:numel(sublist)
    
    substr = sprintf('S%02d',sublist(ss));
    
    behavior_dir = fullfile(exp_path,'DataBehavior',substr);
    output_dir = fullfile(exp_path,'AnalyzeSilhouetteLoc',substr,'EVs');
    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end
    
    sessions = dir([behavior_dir, '/Session*']);
    cd(output_dir)

    for si = 1:numel(sessions)
        session = str2double(sessions(si).name(end));
        loc_files = dir(fullfile(behavior_dir,sessions(si).name, '*Silhouette*'));
        if numel(loc_files)==0
            continue
        end
        assert(numel(loc_files)==1)
        load(fullfile(loc_files(1).folder, loc_files(1).name));
        
        fprintf('found %d runs of silhouette localizer for %s, session %d\n',length(TheData), substr, session);
        
        
        for run = 1:length(TheData)
            
            fprintf('[black white] of checkerboard is [%d %d]\n',TheData(run).p.black,TheData(run).p.white);
            % 12.8s (16TRs*.8s TR) need to be subtracted as they're not in 
            % the volumes (they were calibration or whatever TRs)
            not_recorded = 12.8;     
           
           
            % To make output files for the conditions (one condition, stim on/off). 
            % Each event (row) is described by three numbers: the onset in 
            % seconds (first column), the duration in seconds (second 
            % column), and the value of the input during that period.
            % very last event is the fixation that begins immediately after
            % offset of the last stimulus, has a duration of last iti+5
            trial_starts = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime;
            dur = TheData(run).t.StimTimeTotal;
            onset = [0; trial_starts-not_recorded; (trial_starts(end)-not_recorded)+dur];
            duration = [TheData(run).t.StartFixation-not_recorded; repmat(dur, TheData(run).p.nTrials,1); TheData(run).t.ITI(end)+TheData(run).t.EndFixation];
                
    
            input = [0; ones(TheData(run).p.nTrials,1); 0]; 

            stim_text = [onset, duration, (input==1)];

            % now write the text files
            formatspec = '%.1f % .1f % d\n';

            filename = sprintf('%s_session%d_run%d_stim.txt',substr,session,run);
            h = fopen(filename,'w');
            fprintf(h,formatspec,stim_text');
            fclose(h);


        end %run
        
    end %session
    
end %subject

cd(my_dir)






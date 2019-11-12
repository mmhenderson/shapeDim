%%%%%%%%%%%%%%% FIRST SCRIPT TO RUN FOR LOCALIZER ANALYSES %%%%%%%%%%%%%%%%
% get timing for each localizer run, and save out to text files that can be
% read by FEAT (one file per explanatory viariable per run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


subs = {'01'};
my_dir = pwd;

root = '/usr/local/serenceslab/maggie/shapeDim/Pilot1/';
locfolder = 'AnalyzeCatLocalizer';

for n = 1:numel(subs)
    
    fid = fopen([root 'DataPreproc/S' subs{n} '/runs.list'],'r');
    runstrs = [];
    line = fgetl(fid);
    while ischar(line) && ~isempty(line)
        if contains(line, 'floc')
            runstrs = [runstrs; line(1:6)];
        end
        line = fgetl(fid);
    end
    sess2print = runstrs(:,2);
    runs2print = runstrs(:,6);
    fclose(fid);
    
    if strcmp(subs{n},'01')
        % we lost data for run 1, so use run 2 only
        sess2print = sess2print(2);
        runs2print = runs2print(2);
    end
    
    behavior_dir = [root 'DataBehavior/S',...
        char(subs(n))];
    output_dir = [root locfolder '/S',...
        char(subs(n)), '/EVs'];
    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end
    
    loc_files = dir([behavior_dir, '/*/*fLoc*_alldat.mat']);
    cd(output_dir)

    unvals = str2num(unique(sess2print))';
    for session = unvals
        
        loc_files = dir([behavior_dir, '/Session' sess2print(session) '/*fLoc*alldat.mat']);
        
        fprintf('session %d: found %d category localizer files\n',session, length(loc_files));
            
        for run = 1:length(loc_files)
            
            if ~contains(loc_files(run).name, ['run' num2str(run)])
                error('run number of files doesn''t match')
            end
            
            load([loc_files(run).folder, '/', loc_files(run).name])
         
            % catDirs stores the string labels for each image category-
            % 'scrambled'    'scrambled'    'body'    'limb'    'adult'    'child'    'corridor'    'house'    'car'    'instrument'
            % note that these are in pairs - body and limb are both body
            % parts, etc. To simplify the analysis we will group by pairing
            % now. 
            
            labels = {'scrambled','body','face','place','object','blank'};
            
            % there are this many conditions, plus a blank cond
            numConds=length(labels);
            
            % 12.8s (16TRs*.8s TR) need to be subtracted as they're not in 
            % the volumes (they were calibration or whatever TRs)
            not_recorded = 12.8;     
            
            
            % find indices of all trials where a new block began (whether
            % it's blank, or a category condition).
            trialinds_startblock = [1, find(diff(theSubject.trials.cond)~=0)+1];
            trialinds_stopblock = [find(diff(theSubject.trials.cond)~=0), length(theSubject.trials.cond)];
            
            ntrials_eachblock = diff([trialinds_startblock',trialinds_stopblock'],[],2)+1;
            
            % mark the onset, duration, and condition of each block
            onset = theSubject.trialOnset(trialinds_startblock)' - theSubject.triggerOnset - not_recorded;
            duration = ntrials_eachblock*viewTime;
            
            % here there are 10 total condition labels, but there are
            % actually only 5 categories (2 conditions per category). Make
            % that conversion here:
            cond_10 = theSubject.trials.cond(trialinds_startblock)';
           
            cond = cond_10;
            cond(cond_10==1 | cond_10==2) = 1;
            cond(cond_10==3 | cond_10==4) = 2;
            cond(cond_10==5 | cond_10==6) = 3;
            cond(cond_10==7 | cond_10==8) = 4;
            cond(cond_10==9 | cond_10==10) = 5;
            % blank is category 6
            cond(cond_10==0) = 6;
            
             % now write the text files
            formatspec = '%.1f % .1f % d\n';

            
            for ii=1:numConds
                
                text = [onset duration (cond==ii)];
                
                filename = ['S', char(subs(n)), '_session', num2str(session), '_run', num2str(run) '_' labels{ii} '.txt'];
                h = fopen(filename,'w');
                fprintf(h,formatspec,text');
                fclose(h);
                
            end

        end %run
        
    end %session
    
end %subject

cd(my_dir)






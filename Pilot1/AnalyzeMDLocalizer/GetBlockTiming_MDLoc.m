%%%%%%%%%%%%%%% FIRST SCRIPT TO RUN FOR LOCALIZER ANALYSES %%%%%%%%%%%%%%%%
% get timing for each localizer run, and save out to text files that can be
% read by FEAT (one file per explanatory viariable per run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


subs = {'01'};
my_dir = pwd;

root = '/usr/local/serenceslab/maggie/shapeDim/Pilot1/';
locfolder = 'AnalyzeMDLocalizer';

for n = 1:numel(subs)
    
    fid = fopen([root 'DataPreproc/S' subs{n} '/runs.list'],'r');
    runstrs = [];
    line = fgetl(fid);
    while ischar(line) && ~isempty(line)
        if contains(line, 'mdloc')
            runstrs = [runstrs; line(1:6)];
        end
        line = fgetl(fid);
    end
    sess2print = runstrs(:,2);
    runs2print = runstrs(:,6);
    fclose(fid);
    
    behavior_dir = [root 'DataBehavior/S',...
        char(subs(n))];
    output_dir = [root locfolder '/S',...
        char(subs(n)), '/EVs'];
    if ~exist(output_dir,'dir')
        mkdir(output_dir);
    end
    
%     loc_files = dir([behavior_dir, '/*/*_alldat.mat']);
    cd(output_dir)

    unvals = str2double(unique(sess2print))';
    for session = 1:length(unvals)
        
        loc_files = dir([behavior_dir, '/Session' num2str(unvals(session)) '/*MDLocRun*.mat']);
        
        fprintf('session %d: found %d MD localizer files\n',unvals(session), length(loc_files));
            
        for run = 1:length(loc_files)
            
           
            if ~contains(loc_files(run).name, ['Run' num2str(run)])
                error('run number of files doesn''t match')
            end
            
            load([loc_files(run).folder, '/', loc_files(run).name])
            
            labels = {'Easy','Hard'};
            
            % there are this many conditions, plus a blank cond
            numConds=length(labels);
            
            % 12.8s (16TRs*.8s TR) need to be subtracted as they're not in 
            % the volumes (they were calibration or whatever TRs)
            not_recorded = 12.8;     
            
            %start time of the exp, after 16 frames trimmed
            startExpTrunc=t.TTL_time+not_recorded;
            
            onset = [];
            duration = [];
            cond = [];
            
            for ii=1:numConds

                % in the third column of d.data - 0=blank, 2=easy, 3=hard
                diffs=diff(d.data(:,3)==ii+1);
                %indices of the first trials in this block type
                indsstart=find(diffs==1)+1;
                %indices in the last trials of this block type
                indsstop=find(diffs==-1);

                theseonsets = t.flip_times(indsstart,1) - startExpTrunc;
                thesedurations = t.flip_times(indsstop,7) - t.flip_times(indsstart,1);

                onset = [onset; theseonsets];
                duration = [duration; thesedurations];
                cond = [cond; repmat(ii, length(theseonsets), 1)];
                
            end
            
            [onset,sortorder] = sort(onset, 'ascend');
            duration = duration(sortorder);
            cond = cond(sortorder);
           
             % now write the text files
            formatspec = '%.1f % .1f % d\n';

            
            for ii=1:numConds
                
                text = [onset duration (cond==ii)];
                
                filename = ['S', char(subs(n)), '_session', num2str(unvals(session)), '_run', num2str(run) '_' labels{ii} '.txt'];
                h = fopen(filename,'w');
                fprintf(h,formatspec,text');
                fclose(h);
                
            end

        end %run
        
    end %session
    
end %subject

cd(my_dir)






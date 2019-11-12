%%%%%%%%%%%%%% SECOND SCRIPT TO RUN FOR LOCALIZER ANALYSES %%%%%%%%%%%%%%%%
% To run the GLM this script will first make design .fsf files for each run,
% and next run the localizer analysis (GLM) run for run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


subs = {'01'};
my_dir = pwd;
project_dir = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';
locfolder = 'AnalyzeCatLocalizer';

spatial_smoothing = false; % do or don't do 2mm spatial smoothing
higher_level_FE = true; % do higher level Fixed Effects (if false: does Mixed Effects i.e. FLAME 1).

for n = 1:numel(subs)
    
    %%%% input %%%%
    % where are the Explanatory Variable text files at?
    EV_dir = [project_dir, locfolder '/S', char(subs(n)), '/EVs'];
    % where are my 4D input localizer nifti's at (make sure a runs.list is saved in here, too)?
    nifti_dir = [project_dir, '/DataPreproc/S', char(subs(n))];
    niftiID = 'REG_MC_DET';
    runnrs = [];
    runsid = fopen([nifti_dir,'/runs.list']);
    line = fgetl(runsid);
    while ischar(line) && ~isempty(line)
        if contains(line,'floc')
            runnrs = [runnrs; line(1:6)];
        end
        line = fgetl(runsid);
    end
    fclose(runsid);
    
    if strcmp(subs{n},'01')
        % we lost data for run 1, so use run 2 only
        runnrs = runnrs(2,:);
    end
    
    %%%% output %%%%
    % where will the .fsf files go
    fsf_dir = [project_dir, locfolder '/S', char(subs(n)), '/fsfs'];
    if ~exist(fsf_dir,'dir')
        unix(['mkdir -p ', fsf_dir]);
    end
    % where will the .feat directories go
    feat_dir = [project_dir, locfolder '/S', char(subs(n)), '/feats'];
    if ~exist(feat_dir,'dir')
        unix(['mkdir -p ', feat_dir]);
    end
    
    
    for run = 1:size(runnrs,1)
        
        % run id string
        runidstr = ['S', char(subs(n)), '_', char(runnrs(run,1:2)), '_', char(runnrs(run,4:6))];
        runIDs{run} = runidstr;
        
        % get id string to input 4D functional volume
        idstr_input_4D = dir([nifti_dir, '/' char(runnrs(run,1:2)) '_', '*', niftiID, '*', char(runnrs(run,4:6)), '*.nii.gz']);
        idstr_input_4D = [idstr_input_4D.folder, '/', idstr_input_4D.name(1:end-7)];
        
        % get id string to custom EV files
        runs_this_session = runnrs(runnrs(:,2)==runnrs(run,2),4:6);
        for numruns = 1:size(runs_this_session,1)
            if strcmp(runs_this_session(numruns,:), runnrs(run,4:6))
                loc_run = numruns; %find # nth localizer run in this session 
            end
        end
        
        labels = {'scrambled','body','face','place','object','blank'};
        
        EV_files = [];
        for ii=1:length(labels)
            EV_files = [EV_files, dir([EV_dir, '/*session', char(runnrs(run,2)),...
                '*run', num2str(loc_run), '_' labels{ii} '.txt'])];
            
        end
        
        %         idstr_EV1 = [EV_dir, '/', EV_files(1).name]; % left
%         idstr_EV2 = [EV_dir, '/', EV_files(2).name]; % right
        
        % create name & place for output folder (i.e. for .feat directories
        % to go into)
        idstr_output = [feat_dir, '/', runidstr];
        
        %%%% make .fsf files
        % read my .fsf template and write out a run specific .fsf
        fin = fopen([project_dir, locfolder '/design_template.fsf'],'r');
        fout = fopen([fsf_dir, '/', runidstr, '_design.fsf'],'w+');
        while 1
            line = fgetl(fin);
            % replace output directory
            if contains(line,'set fmri(outputdir)')
                line = ['set fmri(outputdir) "', idstr_output, '"'];
                fprintf([line, '\n'])
            end
            % replace input data
            if contains(line,'set feat_files(1)')
                line = ['set feat_files(1) "', idstr_input_4D, '"'];
                fprintf([line, '\n'])
            end
             
            
            for ii=1:length(EV_files)
            % replace EVs
                if contains(line,['set fmri(custom' num2str(ii) ')'])
                    idstr = [EV_dir, '/', EV_files(ii).name];
                    line = ['set fmri(custom' num2str(ii) ') "', idstr, '"'];
                    fprintf([line, '\n'])
                end
                if contains(line,['set fmri(evtitle' num2str(ii) ')'])                   
                    line = ['set fmri(evtitle' num2str(ii) ') "', labels{ii} , '"'];
                    fprintf([line, '\n'])
                end
            end
            
            % replace spatial smoothing (default is on) if unwanted
            if ~spatial_smoothing
                if contains(line,'set fmri(smooth) 2.0')
                    line = 'set fmri(smooth) 0';
                    fprintf([line, '\n'])
                end
            end
            % print line to new .fsf
            fprintf(fout,[line, '\n']);
            % find the last line and break this loop
            if contains(line,'set fmri(overwrite_yn)')
                break
            end
        end
        fclose(fin);
        fclose(fout);
        
        %%%% run the GLM on this run
        runglm = true;
        if exist([idstr_output '.feat'],'dir')
            fprintf(['A .feat directory for run ',  char(runnrs(run,4:6)), ', session ',...
                char(runnrs(run,2)), ', already exists (with ',...
                num2str(size(dir([idstr_output '.feat']),1)) ,'files in it)...\n']);
            rsp = input('Want to overwrite? (y/n) --> ', 's');
            if strcmp(rsp, 'n')
                runglm = false;
            end
        end
        if runglm
            fprintf(['Running GLM on run ',  char(runnrs(run,4:6)), ' from session ', char(runnrs(run,2)), '...\n'])
            fprintf('...this will take a while...\n')
            unixstr = ['feat ', fsf_dir, '/', runidstr, '_design.fsf'];
            [~,s] = unix(unixstr,'-echo');
            if s == 1 %a non-zero value returned by unix means the command has failed
                fprintf('\n\nERROR: Failed to run feat sucessfully. \n');
            end
        end
        
        %%%% do workaround for future higher-level analysis
        %create 'reg' folder in this runs .feat folder
        doworkaround = false;
        if ~exist([idstr_output, '.feat/reg'],'dir')
            unix(['mkdir ', idstr_output, '.feat/reg']);
            doworkaround = true;
        end
        if doworkaround
            %add example_func2standard.mat and standard2example_func.nii.gz to this folder
            unix(['scp ', project_dir, locfolder '/*example_func*.mat ', idstr_output, '.feat/reg/']);
            %copy example_func.nii.gz from .feat & rename it standard.nii.gz
            unix(['scp ', idstr_output, '.feat/example_func.nii.gz ', idstr_output, '.feat/reg/standard.nii.gz']);
        end
    end %run
    
    
    %% Run higher-level GLM (across all runs)
    %%%% make .fsf file
    % read my .fsf template and write out a run specific .fsf
    hfin = fopen([project_dir, locfolder '/high_level_design_template.fsf'],'r');
    hfout = fopen([fsf_dir, '/S', char(subs(n)), '_high_level_design.fsf'],'w+');
    while 1
        line = fgetl(hfin);
        % replace output directory
        if contains(line,'set fmri(outputdir)')
            line = ['set fmri(outputdir) "', feat_dir, '/AllSessionsFE"'];
            fprintf([line, '\n'])
        end
        % set stats effects (FE or ME)
        if contains(line,'set fmri(mixed_yn)')
            if ~higher_level_FE % default if FE (i.e. 3)
                line = 'set fmri(mixed_yn) 2';
                fprintf([line, '\n'])
            end
        end
        %replace input data
        for run = 1:size(runnrs,1)
            lookfor = ['set feat_files(', num2str(run), ')'];
            if contains(line,lookfor)
                line = [lookfor, ' "', feat_dir, '/', char(runIDs(run)), '.feat"'];
                fprintf([line, '\n'])
            end
        end
        % print this line to new .fsf
        fprintf(hfout,[line, '\n']);
        % find the last line and break this loop
        if contains(line,'set fmri(overwrite_yn)')
            break
        end
    end
    fclose(hfin);
    fclose(hfout);
    %%%% run the GLM on this run
    fprintf('Running high-level GLM across all runs...\n')
    unixstr = ['feat ', fsf_dir, '/S', char(subs(n)), '_high_level_design.fsf'];
    [~,s] = unix(unixstr,'-echo');
    if s == 1 %a non-zero value returned by unix means the command has failed
        fprintf('\n\nERROR: Failed to run feat sucessfully. \n');
    end
    
end %subject







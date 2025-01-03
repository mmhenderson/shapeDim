%%%%%%%%%%%%%% SECOND SCRIPT TO RUN FOR LOCALIZER ANALYSES %%%%%%%%%%%%%%%%
% To run the GLM this script will first make design .fsf files for each run,
% and next run the localizer analysis (GLM) run for run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all


sublist = [8,9,10];
my_dir = pwd;
% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 1;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));

subfolder = 'AnalyzeSilhouetteLoc';
nStims = 1;    % how many stim types (predictors)
% how many TRs in this kind of run, after trimming?
% nTRs = 328 - 16;

spatial_smoothing = false; % do or don't do 2mm spatial smoothing
higher_level_FE = true; % do higher level Fixed Effects (if false: does Mixed Effects i.e. FLAME 1).

% usually set these both to 1, but maybe you got stuck somewhere in the
% middle and need to do HL only.
doLL = 1;
doHL = 1;
% note, the HL analysis only works from this script if your number of localizer
% runs exactly matches what is in the high_level_design_template.fsf file.
% It will throw an error if not - in that case you will need to either
% create a new template yourself (can edit it as a text file or in the GUI)
% and re-run this script, or you can bypass this script entirely and run
% the analysis straight from the GUI.

for ss = 1:numel(sublist)
    
    if sublist(ss)<8
        % GE setup
        nTRs = 328 - 16;
        tr_dur = 0.8;
    else
        % Prisma setup
        nTRs = 202;
        tr_dur = 1.3;
    end
    substr = sprintf('S%02d',sublist(ss));
    fprintf('%s, expect %d TRs\n', substr, nTRs)
    
    % where are the EV files?
    EV_dir = fullfile(exp_path, subfolder, substr, 'EVs');
    
    % where will FSFs go?
    fsf_dir = fullfile(exp_path, subfolder, substr, 'fsfs');
    if ~exist(fsf_dir,'dir')
        unix(['mkdir -p ', fsf_dir]);
    end
    
    % where will the .feat directories go
    feat_dir = fullfile(exp_path, subfolder, substr, 'feats');
    if ~exist(feat_dir,'dir')
        unix(['mkdir -p ', feat_dir]);
    end
    
    nifti_dir = fullfile(exp_path, 'DataPreproc',substr);
   
    niftiID = 'REG_MC_DET';
    % where are my 4D input localizer nifti's at (make sure a runs.list is saved in here, too)?
    runnrs = [];
    runsid = fopen([nifti_dir,'/runs.list']);
    line = fgetl(runsid);
    while ischar(line) && ~isempty(line)
        if contains(line,'silh')
            runnrs = [runnrs; line(1:6)];
        end
        line = fgetl(runsid);
    end
    fclose(runsid);
    
    runIDs = [];
    for run=1:size(runnrs,1)
        % run id string
        runidstr = [substr, '_', char(runnrs(run,1:2)), '_', char(runnrs(run,4:6))];
        runIDs{run} = runidstr;
    end
    
    if doLL
    
        for run = 1:size(runnrs,1)

            runidstr = runIDs{run};

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
            
           
            EV_files = [];
            labels = {'stim'};
            
            EV_files = [EV_files, dir([EV_dir, '/*session', char(runnrs(run,2)),...
                '*run', num2str(loc_run), '_stim.txt'])];
             

            % create name & place for output folder (i.e. for .feat directories
            % to go into)
            idstr_output = [feat_dir, '/', runidstr];

            %%%% make .fsf files
            % read my .fsf template and write out a run specific .fsf
            fin = fopen(fullfile(exp_path,subfolder, 'design_template.fsf'),'r');
            fout = fopen(fullfile(fsf_dir, [runidstr, '_design.fsf']), 'w+');
            
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
                % replace EVs
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
                
                if contains(line, 'set fmri(tr)')
                    
                    line = ['set fmri(tr) ' num2str(tr_dur)];
                    
                    fprintf([line '\n']);
                end
                
                if contains(line, 'set fmri(npts)')
                    
                    line = ['set fmri(npts) ' num2str(nTRs)];
                    
                    fprintf([line '\n']);
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
                f2copy = fullfile(exp_path,subfolder,'*example_func*.mat');
                unix(['scp ', f2copy, ' ', idstr_output, '.feat/reg/']);
                %copy example_func.nii.gz from .feat & rename it standard.nii.gz
                unix(['scp ', idstr_output, '.feat/example_func.nii.gz ', idstr_output, '.feat/reg/standard.nii.gz']);
            end
        end %run

    end
    
    %% Run higher-level GLM (across all runs)
    
    if doHL
       

        %%%% make .fsf file
        % read my .fsf template and write out a run specific .fsf
        hfin = fopen(fullfile(exp_path,subfolder,sprintf('high_level_design_template_%d.fsf',size(runnrs,1))), 'r');
      
        % figure out which template file matches my number of runs.
        hfout = fopen(fullfile(fsf_dir, [substr, '_high_level_design.fsf']), 'w+');
      
        if hfin==-1
            error('No template fsf file to match the number of runs you did (%d). Use the FEAT GUI!', size(runnrs,1))
        end
            
        % a flag to verify that we have the correct number of runs...
        num_runs_checked = 0;
        while 1
            line = fgetl(hfin);
            if line==-1
                break
            end
            if contains(line, 'set fmri(npts)')
                try
                    num_runs = str2double(line(16:17));
                catch
                    num_runs = str2double(line(16));
                end
                if num_runs~=size(runnrs,1)
                    error('oops, check npts in your FSF file')
                else
                    num_runs_checked = 1;
                end
            end
            
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
        
        if num_runs_checked==0
            error('could not confirm that your fsf file has the correct number of inputs. exiting.')
        end
        %%%% run the GLM on this run
        fprintf('Running high-level GLM across all runs...\n')
        unixstr = ['feat ' fullfile(fsf_dir, [substr, '_high_level_design.fsf'])];
        
        [~,s] = unix(unixstr,'-echo');
        if s == 1 %a non-zero value returned by unix means the command has failed
            fprintf('\n\nERROR: Failed to run feat sucessfully. \n');
        end

    end
end %subject







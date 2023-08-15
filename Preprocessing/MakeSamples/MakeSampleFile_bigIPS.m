function [] = MakeSampleFile_bigIPS(sub2do, debug)

%% Make Sample File for shapedim
% This will be a big .mat file that concatenates data from all runs of 
% each task, each having the format [nTRsTotal x nVoxels]. 
% All preprocessing of task data and analysis of the spatial localizer
% (used for thresholding) are run before this. 

%%

close all

if nargin<2
    debug=false;
end

if nargin<1
    sub2do = [1,2,3,4,5,6,7]
end

subinit_big = {'CP','BX','BR','CA','CG','CT','CU'};
subnum_big = [1,2,3,4,5,6,7];

% find my root directory - up a few dirs from where i am now
mypath = pwd;
filesepinds = find(mypath==filesep);
nDirsUp = 2;
exp_path = mypath(1:filesepinds(end-nDirsUp+1));
hemis = {'lh', 'rh'};   
addpath('/usr/local/freesurfer6/matlab/')   % this folder has the load_nifti function we need

% which areas do we want to load voxels from?
ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

% the actual filenames for the ROIs will have underscores instead of spaces
ROI_fns = ROI_names;
for aa=1:numel(ROI_fns)
    ROI_fns{aa}(ROI_fns{aa}==' ') = '_';
end
nVOIs = numel(ROI_names);

% set paths 
out_path = fullfile(exp_path, 'Samples');
beh_path = fullfile(exp_path, 'DataBehavior');

for ss = sub2do
    
    % get subject information
    subinit = subinit_big{ss};
    substr = sprintf('S%02d',subnum_big(ss));
    
    % set subject preproc data path
    func_path = fullfile(exp_path, 'DataPreproc',substr);  
    % this is where the t-map from the localizer lives
    feat_path = fullfile(exp_path, 'AnalyzeSilhouetteLoc',substr,'feats','AllSessionsFE.gfeat','cope1.feat');
   
    % the file I will save out.
    fn2save = fullfile(out_path,sprintf('SampleFile_bigIPS_%s.mat',substr));
   
    %% Localizer stats
     %Get the t-stat maps that I got from doing the GLM on my localizer(s)
    fprintf('Loading localizer thresholding map from %s\n', fullfile(feat_path, 'stats','tstat1.nii.gz'));
    nifti = load_nifti(fullfile(feat_path, 'stats','tstat1.nii.gz'));
    vmpVals = nifti.vol; %Store t-values for each voxel index that is in my visual areas in 'vmpVals'
    
    cd(fullfile(feat_path,'stats'));
    % get FE dof
    [~, my_FE_dof] = unix('fslstats tdof_t1 -M'); % gives the dof I need for fixed effects (which is the test I ran), output is a string
    
    % make log p-stats map
    unix(['ttologp -logpout logp1 varcope1 cope1 ', num2str(str2num(my_FE_dof))]);
    
    % convert from log to p-value map
    unix(['fslmaths logp1 -exp p1']);
    
    % do FDR on p-value map and get probability threshold
    [~,prob_thresh] = unix('fdr -i p1 -q 0.05');

    % go from my p-threshold back into my t-stat
    my_t = abs(tinv(str2num(prob_thresh(28:end)),str2num(my_FE_dof)));
    
    tThresh = my_t;
    vox_above_thresh = find(vmpVals>tThresh);
    
    %% Read VOIs
    % Read the retino VOI files, which consist of a 3D volume with 1s marking
    % the location of the ROI. We will store all the indices of each ROI, and
    % then make sure none overlap.
    % Also loading all my motor ROIs here. All motor ROIs will be in the VOIs
    % folder, with the suffix DIGITLOC to indicate how they were generated.
    % Later we will threshold the retinotopic ROIs, but not the digit ROIs,
    % with our orientation localizer - but that does not happen in this script.
    
    fprintf('Loading all ROIs\n')
    
    % this will be a structure array [2 x nROIs] listing all my rois and their properties.
    ROIs = struct('voxel_inds',[], 'name',' ','is_visual',[],'is_motor',[],'is_md',[]);
    
    % loop over hemispheres
    for hh = 1:length(hemis)
        
        % loop over VOIs
        for vv = 1:nVOIs
            
            % looking for VOI files with this name...
            if vv<4
                % early visual area, has a dorsal and ventral component
                file1 = dir(fullfile(exp_path,'VOIs',substr,sprintf('%s_%sd.nii.gz',hemis{hh},ROI_fns{vv})));
                file2 = dir(fullfile(exp_path,'VOIs',substr,sprintf('%s_%sv.nii.gz',hemis{hh},ROI_fns{vv})));
                VOIfiles = [file1;file2];
                assert(numel(VOIfiles)==2);
            else
                % only has one part.
                file1 = dir(fullfile(exp_path,'VOIs',substr,sprintf('%s_%s.nii.gz',hemis{hh},ROI_fns{vv})));
                VOIfiles = file1;
                assert(numel(VOIfiles)==1);
            end
           
            % if it's missing, print a warning here but keep going. 
            if numel(VOIfiles)==0
                fprintf('missing %s %s for %s\n', hemis{hh}, ROI_names{vv}, substr);
                continue
            end
           
            % check if this is a motor ROI
            ROIs(hh,vv).is_motor = false;  
            ROIs(hh,vv).is_md = false;  
            ROIs(hh,vv).is_visual = true;
            
            % now load the niftis, loop over individual parts of the ROI 
            % (dorsal/ventral if applicable)  
            voxel_inds_this_ROI = [];
            for ff = 1:length(VOIfiles)

                nifti = load_nifti(fullfile(VOIfiles(ff).folder,VOIfiles(ff).name));
            
                % extracting the indices from the binary mask with "find". now each
                % voxel gets a unique number that will stay with it throughout
                % preprocessing.
                voxel_inds_this_ROI = [voxel_inds_this_ROI; find(nifti.vol)];

            end
            % if there are multiple parts, get rid of duplicates
            voxel_inds_this_ROI =unique(voxel_inds_this_ROI);
            
            % of these indices, which ones are above threshold from
            % localizer mask?
            do_thresh = ~contains(ROI_names{vv}, 'IPS');
            if do_thresh
                inds2use = intersect(voxel_inds_this_ROI, vox_above_thresh);
                fprintf('%s %s %s: %.2f percent vox are above spatial localizer threshold\n',substr,hemis{hh},ROI_names{vv},numel(inds2use)/numel(voxel_inds_this_ROI)*100);
            else
                inds2use = voxel_inds_this_ROI;
                fprintf('%s %s %s: using all voxels\n',substr,hemis{hh},ROI_names{vv});
            end
           
            % add information to the structure
            ROIs(hh,vv).voxel_inds = inds2use';
            ROIs(hh,vv).name = sprintf('%s_%s',hemis{hh},ROI_names{vv});
           
        end
    end
  
    
    %% Correct overlap
    % Prune voxels which are shared between pairs of ROIs (interleaved):
    for hh = 1:length(hemis)
        for vv1 = 1:nVOIs            
            %for the n visual areas in this hemisphere, compare that visual area against all other n-1 areas
            for vv2 = vv1+1:nVOIs
                
                % find the overlap
                overlapvox = intersect(ROIs(hh,vv1).voxel_inds, ROIs(hh,vv2).voxel_inds);
                if ~isempty(overlapvox)
                    fprintf('detected %d voxels of overlap between %s and %s\n',numel(overlapvox),ROIs(hh,vv1).name,ROIs(hh,vv2).name);
                    
                    if (ROIs(hh,vv1).is_visual + ROIs(hh,vv2).is_visual) == 1    
                       
                        % if one is motor/MD and the other is visual - keep all
                        % voxels in the visual ROI (should only come up for
                        % occipital motor ROIs that we may find, and these
                        % probably won't get used anyway).
                        if ROIs(hh,vv2).is_visual
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv2).name);
                            roi1_deletevox = overlapvox;
                            roi2_deletevox = [];
                        else
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv1).name);
                            roi1_deletevox = [];
                            roi2_deletevox = overlapvox;
                        end
                    elseif (ROIs(hh,vv1).is_md + ROIs(hh,vv2).is_md) == 1  
                        % one is MD, the other is motor - put all voxels in
                        % the motor ROI
                        if ROIs(hh,vv2).is_motor
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv2).name);
                            roi1_deletevox = overlapvox;
                            roi2_deletevox = [];
                        else
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv1).name);
                            roi1_deletevox = [];
                            roi2_deletevox = overlapvox;
                        end
                    elseif (numel(ROIs(hh,vv1).voxel_inds)<5 && numel(ROIs(hh,vv2).voxel_inds)>5) ||...
                            (numel(ROIs(hh,vv2).voxel_inds)<5 && numel(ROIs(hh,vv1).voxel_inds)>5)
                        % if one is very small and other is not, let the
                        % little one keep all its voxels
                        if numel(ROIs(hh,vv2).voxel_inds)<5
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv2).name);
                            roi1_deletevox = overlapvox;
                            roi2_deletevox = [];
                        elseif numel(ROIs(hh,vv1).voxel_inds)<5
                            fprintf('    putting all overlap voxels in %s\n',ROIs(hh,vv1).name);
                            roi1_deletevox = [];
                            roi2_deletevox = overlapvox;
                        end
                    else
                        fprintf('    splitting up evenly\n');                 
                        % otherwise just split them up evenly.
                        roi1_deletevox = overlapvox(1:2:end); %uneven voxels will deleted from roi1
                        roi2_deletevox = overlapvox(2:2:end); %even voxels will be deleted from roi2
                    end
                    ROIs(hh,vv1).voxel_inds(ismember(ROIs(hh,vv1).voxel_inds, roi1_deletevox)) = [];
                    ROIs(hh,vv2).voxel_inds(ismember(ROIs(hh,vv2).voxel_inds, roi2_deletevox)) = [];
                end
            end
        end
    end
    
    % Prune voxels which are shared between hemispheres. For most areas,
    % will just split half into each hemisphere since we'll put them back
    % together later. For motor areas, we care about left/right hemisphere
    % differences, so ditch any voxels that are in both because they're not
    % fully in either hemisphere. 
    overlapvox = intersect([ROIs(1,~[ROIs(1,:).is_motor]).voxel_inds], [ROIs(2,~[ROIs(1,:).is_motor]).voxel_inds]);
    hemi1_deletevox = overlapvox(1:2:end);  % each voxel lives in only one hemisphere now
    hemi2_deletevox = overlapvox(2:2:end);
    fprintf('correcting hemisphere overlap for visual regions: found %d voxels\n',numel(overlapvox))
    
    overlapvox = intersect([ROIs(1,[ROIs(1,:).is_motor]).voxel_inds], [ROIs(2,[ROIs(1,:).is_motor]).voxel_inds]);
    hemi1_deletevox = [hemi1_deletevox, overlapvox];    % all these voxels get ditched
    hemi2_deletevox = [hemi2_deletevox, overlapvox];
    fprintf('correcting hemisphere overlap for motor regions: found %d voxels\n',numel(overlapvox))
    % delete these voxels from whichever ROI they belong to.
    for vv = 1:nVOIs
        
        todelete = intersect(ROIs(1,vv).voxel_inds, hemi1_deletevox);
        if ~isempty(todelete)
            ROIs(1,vv).voxel_inds(ismember(ROIs(1,vv).voxel_inds, todelete)) = [];
        end
        
        todelete = intersect(ROIs(2,vv).voxel_inds, hemi2_deletevox);
        if ~isempty(todelete)
            ROIs(2,vv).voxel_inds(ismember(ROIs(2,vv).voxel_inds, todelete)) = [];
        end
    end
    
    % finally, unwrap all voxels that landed in any of these ROIs.
    all_vox_concat = [ROIs.voxel_inds];
    all_vox_concat = unique(all_vox_concat);
    
    
    %% Get Sample Timecourse
    % This will be a matrix of (timepoints x voxels) big, so for each TR (i.e.
    % length(EventLabels)) by each visually driven voxel (i.e. length(vInd))
    
    fprintf('Loading niftis and extracting time courses\n')
    
    TRditched = 16;
    
    nTRs_main = 327- TRditched;
    nTRs_rep = 329 - TRditched;
   
    RunsListMain = []; RunsListRep = [];
    main_sess_list = [];
    main_run_in_sess_list = [];
    main_rc=0;
    main_se=1;
    fid = fopen(fullfile(func_path,'runs.list'));
    line = fgetl(fid);

    while ischar(line) && ~isempty(line)
        if strcmp(line(8:end),'main')
            RunsListMain = [RunsListMain; line(1:6)];
            se = str2double(line(1:2));
            if se~=main_se
                main_rc = 1;
                main_se = se;
            else
                main_rc = main_rc+1;
            end
            main_sess_list = [main_sess_list; main_se];
            main_run_in_sess_list = [main_run_in_sess_list; main_rc];
            
        end
        if strcmp(line(8:end),'rep')
            RunsListRep = [RunsListRep; line(1:6)];
        end
        
        line = fgetl(fid);
        
    end
    fclose(fid);
    
    % these are files that i made manually, listing which runs we need to
    % load for the main task. For some subjects, we're loading a run from
    % session 3 in place of one from session 1, etc. so this is addressed
    % in this file...
    main_run_sequence = csvread(fullfile(out_path,'run_sequences',sprintf('%s_main_run_sequence.csv',substr)),1,0);
    main_sess_num_load = main_run_sequence(:,4);
    main_run_in_sess_load = main_run_sequence(:,5);
%     main_sess_date_load = main_run_sequence(:,6);
    
    RunsListMainOrig = RunsListMain;
    RunsListMain = [];
    for rr = 1:length(main_sess_num_load)
        idx = (main_sess_list==main_sess_num_load(rr)) & ...
            (main_run_in_sess_list==main_run_in_sess_load(rr));
        run_str = RunsListMainOrig(idx,:);
        RunsListMain = [RunsListMain; run_str];
    end
    assert(length(main_sess_num_load)==length(RunsListMain))
    any_change = (length(RunsListMain)~=length(RunsListMainOrig)) || ...
                ~all(all(RunsListMainOrig==RunsListMain));
    if any_change
        assert(strcmp(substr, 'S02') || strcmp(substr, 'S07'))
        fprintf('Loading main task runs for %s:\n', substr)
        disp(RunsListMain)
    end
    
%     continue
    
    if ~debug
        % Make samplesMain
        fprintf('    main task niftis...\n')
        tmp=NaN(nTRs_main,length(all_vox_concat),length(RunsListMain));
        for run = 1:length(RunsListMain) %for each functional run load the corresponding nifti
            niifile = fullfile(func_path,[RunsListMain(run,1:2), '_REG_MC_DET_', RunsListMain(run,4:6), '.nii.gz']);
            fprintf('loading from %s\n',niifile);
            nifti_run = load_nifti(niifile);
            reshaped_volume = reshape(nifti_run.vol, [prod(nifti_run.dim(2:4)) nifti_run.dim(5)])'; %Reshape so that we have one index for each voxel (corresponding to the indices in our VOI file), with voxels as columns
            tmp(:,:,run)=reshaped_volume(:,all_vox_concat);
        end
        samplesMain = [];
        for run = 1:length(RunsListMain)
            samplesMain=[samplesMain;tmp(:,:,run)];
        end
        clear tmp

        % Make samplesRep
        fprintf('    one-back task niftis...\n')
        tmp=NaN(nTRs_rep,length(all_vox_concat),length(RunsListRep));
        for run = 1:length(RunsListRep) %for each functional run load the corresponding nifti
            niifile = fullfile(func_path,[RunsListRep(run,1:2), '_REG_MC_DET_', RunsListRep(run,4:6), '.nii.gz']);
            fprintf('loading from %s\n',niifile);
            nifti_run = load_nifti(niifile);
            reshaped_volume = reshape(nifti_run.vol, [prod(nifti_run.dim(2:4)) nifti_run.dim(5)])'; %Reshape so that we have one index for each voxel (corresponding to the indices in our VOI file), with voxels as columns
            tmp(:,:,run)=reshaped_volume(:,all_vox_concat);
        end
        samplesRep = [];
        for run = 1:length(RunsListRep)
            samplesRep=[samplesRep;tmp(:,:,run)];
        end
        clear tmp
    else
        samplesMain = [];
        samplesRep = [];
    end
    
    %% Save Sample File
    % saving all VOI labels and data here.
   
    if ~exist(out_path, 'dir'), mkdir(out_path); end
    fprintf('saving to %s\n',fn2save);
    save(fn2save, 'ROIs','all_vox_concat',...
        'samplesMain',...
        'samplesRep',...
        '-v7.3');
        
    cd(mypath)

end
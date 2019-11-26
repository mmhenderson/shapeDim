% Define ROIs that show selectivity for hard>easy with MD localizer

% will load the anatomical definitions of a parcellation of the MD network 
% from Fedorenko, Duncan , & Kanwisher, 2013
% downloaded from imaging.mrc-cbu.cam.ac.uk/imaging/MDsystem

% Need to run GLM using feat to create the t-map first.

%%
clear 
close all

subnum = '01';
subinit = 'CG';

% these two contrasts have to correspond to copes 1 and 2!
% contrasts = {'left>right','right>left'};
hemis = {'lh','rh'};
nHemis = length(hemis);

retfolder = '/usr/local/serenceslab/Doreti/';
subretfolder = [retfolder 'ANAT/' subinit '/'];
exp_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';
func_path = [exp_path, 'DataPreproc/S', char(subnum), '/'];

% will need this file later
if ~exist([exp_path 'DataPreproc/S' subnum '/regheadercenterMOD.mat'],'file')
    error('need your registration matrix first!')
end

%% load and transform the TAL parcels
% this is the nifti we got from Fedorenko group, need a copy in this folder
tal_parcels_orig = fullfile(exp_path, 'Processing/MDROI.nii.gz');
tal_brain_orig_mgz = '/mnt/neurocube/local/freesurfer6/average/mni305.cor.mgz';
tal_brain_orig_nii = fullfile(exp_path,'Processing','mni305.cor.nii.gz');
if ~exist(tal_brain_orig_nii,'file')
    err = unix(['mri_convert ' tal_brain_orig_mgz ' ' tal_brain_orig_nii]);
    if err
        error('your mri_convert command failed!')
    end
end
tal_brain_xfm = fullfile(exp_path,'Processing/TAL_brain_xfm.nii.gz');

% this is the registration file that is generated by freesurfer during
% recon-all, maps from the average TAL brain to the subject's anatomical
% scan
tal_reg_file_xfm = fullfile(subretfolder, 'mri','transforms', 'talairach.xfm');
% tal_reg_file = fullfile(exp_path,'Processing','tal2anat_COPIED.dat');

% tal_reg_file = fullfile(exp_path,'Processing','anat2tal.mat');
% 
% 
% if ~exist(tal_reg_file,'file')
%     err = unix(['lta_convert -inmni ' tal_reg_file_xfm ' -outfsl ' tal_reg_file]);
% %     err = unix(['convert_xfm -omat ' tal_reg_file ' ' tal_reg_file_xfm]);
%     if err
%         error('your transformation failed')
%     end
% end
% reverse it
% tal_reg_file_inverse = fullfile(exp_path,'Processing','tal2anat.xfm');
% if ~exist(tal_reg_file_inverse,'file')
%     err = unix(['convert_xfm -omat ' tal_reg_file_inverse ' -inverse ' tal_reg_file]);
%     if err
%         error('your transformation inversion failed')
%     end
% end

% Need a reference volume, this is the anatomical scan
reference_volume_mgz = fullfile(subretfolder, 'mri','orig.mgz');
% make sure we have this in nifti format, if not then convert now
reference_volume_nii = fullfile(subretfolder, 'mri','orig.nii.gz');
if ~exist(reference_volume_nii,'file')
    err = unix(['mri_convert ' reference_volume_mgz ' ' reference_volume_nii]);
    if err
        error('your mri_convert command failed!')
    end
end

%% make the transformation (with a random other transformation file)

% name the output file
xfm_out = fullfile(exp_path,'Processing/MDROI_XFM_TEST.nii.gz');
reg_file = tal_reg_file_xfm;
% reg_file = fullfile(exp_path,'DataPreproc',['S',subnum], 'registration_auto01.mat');
dummy_nii = fullfile(exp_path,'DataPreproc',['S',subnum], 'MCTemplate01.nii.gz');

err =  unix(['flirt -in ', dummy_nii ' -ref ', dummy_nii,...
    ' -applyxfm -init ', reg_file ' -out ', xfm_out]);

%% make the transformation

% name the output file
tal_parcels_xfm = fullfile(exp_path,'Processing/MDROI_XFM.nii.gz');

err =  unix(['flirt -in ', tal_parcels_orig ' -ref ', reference_volume_nii,...
    ' -applyxfm -init ', tal_reg_file_xfm ' -out ', tal_parcels_xfm]);

%% making a manual registration from tal to anat, not perfect

% reg_file_save = fullfile(exp_path,'Processing','tal2anat.dat');
% 
% err = unix(['tkregister2 --s ' subinit ' --targ ' reference_volume_nii ' --mov ' tal_brain_orig_mgz ' --regheader-center --reg ' reg_file_save]);



% % tal_reg_file = fullfile(exp_path,'Processing','tal2anat_COPIED.dat');
% tal_reg_file = fullfile(exp_path,'Processing','tal2anat.dat');
% tal_reg_file = fullfile(subretfolder, 'mri','transforms', 'talairach.xfm');
%%
tal_reg_file = fullfile(exp_path,'Processing','tal2anat_COPIED.dat');

tal_orig = load_nifti(tal_brain_orig_nii);
anat = load_nifti(reference_volume_nii);

parcel_file = load_nifti(tal_parcels_xfm);
parcels = parcel_file.vol;
unvals = unique


%%

err =  unix(['flirt -in ', tal_brain_orig_nii ' -ref ', reference_volume_nii,...
    ' -applyxfm -init ', tal_reg_file ' -out ', tal_brain_xfm ' -v']);
if err
    error('your flirt command failed!')
else
    disp('SUCCESS')
end

%%
tal_reg_file = fullfile(exp_path,'Processing','tal2anat_COPIED.dat');
% tal_reg_file = fullfile(subretfolder, 'mri','transforms', 'talairach.xfm');

err = unix(['tkregister2 --s ' subinit ' --targ ' reference_volume_mgz ' --mov ' tal_brain_orig_mgz ' --reg ' tal_reg_file]);
% err = unix(['tkregister2 --s ' subinit ' --fstal --no-zero-cras']);
%%
file2save = [exp_path 'DataPreproc/' sub '/BrainMask_REG2FUNC.nii.gz'];

fprintf('loading your original brain mask from %s\n',reference_volume_nii)

fprintf('saving your transformed brain mask to %s\n',file2save)

% make the transformation
err =  unix(['flirt -in ', reference_volume_nii ' -ref ', reference_volume,...
    ' -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', file2save]);
if err
    error('your flirt command failed!')
else
    disp('SUCCESS')
end


% and unzip it if you want to
if unzip
    fprintf('unzipping your transformed brain mask...\n')

    err = unix(['gzip -dk ' file2save]);
    if err
        error('your unzip command failed!')
    else
        disp('SUCCESS')
    end
end
%%

% load the annotation file so that we can get a list of all the names (this
% can really be any subject or hemisphere because they're all the same
% names)
% next we will make a volume for every single
% cortical area, then see how many intersect with the face>place localizer. 
[inds,labs,info]= read_annotation('/usr/local/serenceslab/Doreti/ANAT/BX/label/lh.aparc.annot');
allnames = info.struct_names(2:end);
ROIs = allnames;
ROIs_4unix = ROIs;
ROIs_4saving = ROIs;
for rr=1:length(ROIs)
    if contains(ROIs{rr},'&')
        % put an extra slash in here so that unix won't be confused by the & signs
        newname = ROIs{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '\' newname(andind:end)];
        ROIs_4unix{rr} = newname;
        newname = ROIs{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '_and_' newname(andind+1:end)];
        ROIs_4saving{rr} = newname;
    end
end
clear info


% cope 1 is face>place, cope 7 is object>scrambled
cope2use = [3];
copenames = {'mdloc'};

% this first dir is a place to store the un-thresholded VOIs after they are
% made, the second dir is our usual VOI folder
VOIdir1 = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
VOIdir2 = [exp_path, 'VOIs/S' char(subnum) '/'];

minVox = 50;

currentdir = pwd;

% set this to 1 if you want to make all the volumes in functional space.
% takes a while!
doVolumes = 1;

%% first go into the Doreti folder, and make labels based on the anatomy

if doVolumes

    % check to see if labels are already made or not
    makelabs = true;
    if exist([subretfolder, 'label/parcels_all'], 'dir')
        if ~isempty(dir([subretfolder, 'label/parcels_all/*.label'])) % if there are labels in here already
            rsp = input('Labels already exist in this subject''s parcels_all folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makelabs = false;
            end
        end
    else
        mkdir([subretfolder, 'label/parcels_all']);
    end

    if makelabs
        cd(retfolder)
        for hh=1:nHemis
%             unix(['mri_annotation2label --annotation aparc.a2009s --subject ' subinit ' --hemi ' hemis{hh} ' --outdir ' subretfolder 'label/parcels_all/']); 
            unix(['mri_annotation2label --annotation aparc --subject ' subinit ' --hemi ' hemis{hh} ' --outdir ' subretfolder 'label/parcels_all/']); 

        end
        cd(currentdir)
    end
    %% now turn them into VOI masks in functional space

    % here we are making a new folder to hold ROIs that were generated using
    % fsl parcellation scheme. These have to be processed a bit more before
    % they can be put into our normal pipeline. 

    % check to see if masks are already made or not
    makemasks = true;
    if exist([exp_path, 'VOIs_parcel/S', char(subnum)], 'dir')
        if ~isempty(dir([exp_path, 'VOIs_parcel/S', char(subnum), '/*h_*nii.gz'])) % if there are visual masks found
            rsp = input('Visual area masks already exist in this subject''s VOIs_parcel folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makemasks = false;
            end
        end
    else
        mkdir([exp_path, 'VOIs_parcel/S', char(subnum)]);
    end

    % make Visual Area masks
    if makemasks
        delete([exp_path, 'VOIs_parcel/S' char(subnum) '/*h_*nii.gz']); %delete retino content from VOI folder 
    %     ROIs_orig = {'S_precentral-sup-part','S_precentral-inf-part','S_central','S_postcentral',...
    %         'G_precentral','G_postcentral','G_parietal_sup',...
    %         'G&S_subcentral','G_pariet_inf-Supramar','G&S_paracentral','G_front_sup'};

        func_template = [func_path, 'MCTemplateXFM01.nii.gz'];
        regfile = [func_path, 'Func2Anat_auto.dat'];
        outdir = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
        if ~isdir(outdir)
            mkdir(outdir)
        end
        if ~exist(regfile,'file')
            error('registration to anat not made yet! you need to run getVisualROIs.m first')
        end
        % below is in replacement of label2VOI function
        for hemi_idx = 1:length(hemis)
            for roi_idx = 1:length(ROIs_4unix)
                labelfile = [subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', ROIs_4unix{roi_idx}, '.label'];


                if exist([subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', ROIs{roi_idx}, '.label'], 'file')

                    unixcmd = ['mri_label2vol --label ', labelfile, ' --temp ', func_template, ' --subject ', subinit, ' --hemi ', hemis{hemi_idx}, ...
                        ' --o ', outdir, hemis{hemi_idx}, '_', ROIs_4unix{roi_idx}, '.nii.gz --proj frac 0 1 0.1 --fillthresh 0.3 --reg ', regfile];
                    [s,~] = unix(unixcmd, '-echo');
                    if s==1
                        out(['\nERROR: An error occured while processing label ', labelfile, '. Exiting.\n']);
                        status = 1;
                        return
                    end                
                end            
            end %roi_idx   
        end
    end    
end

%% now load the localizer, and mask out the parts we need

for cc=1:length(cope2use)

    stats_path = [exp_path 'AnalyzeMDLocalizer/S' char(subnum) '/feats/AllSessionsFE.gfeat/cope' num2str(cope2use(cc)) '.feat/stats/'];
    cd(stats_path)
%     % load my t-stat map for all voxels
    nifti = load_nifti('tstat1.nii.gz');
    vmpVals = nifti.vol; %Store t-values for each voxel index that is in my visual areas in 'vmpVals'

    clear nifti
    %% find my threshold 
    % this is all copied from Find_t_thresh.m

     % get FE dof
    [~, my_FE_dof] = unix('fslstats tdof_t1 -M'); % gives the dof I need for fixed effects (which is the test I ran), output is a string

    % make log p-stats map
    unix(['ttologp -logpout logp1 varcope1 cope1 ', num2str(str2double(my_FE_dof))]);

    % convert from log to p-value map
    unix(['fslmaths logp1 -exp p1']);

    % do FDR on p-value map and get probability threshold
    [~,prob_thresh] = unix('fdr -i p1 -q 0.05');

    % go from my p-threshold back into my t-stat
    my_t = abs(tinv(str2double(prob_thresh(28:end)),str2double(my_FE_dof)));

    %% define my final mask

    vals2use = vmpVals>my_t;
    
    %% create each VOI file

    for hh=1:length(hemis)

        for vv=1:length(ROIs)

            if exist([subretfolder, '/label/parcels_all/', hemis{hh}, '.', ROIs{vv}, '.label'], 'file')

                fullVOI = [VOIdir1 hemis{hh} '_' ROIs_4unix{vv} '.nii.gz'];
                thisVOI = load_nifti(fullVOI);

                fullVol = thisVOI.vol;

                maskedVol = fullVol & vals2use;
%                 fprintf('%s %s %s: found %d voxels for %s\n',subinit,hemis{hh}, ROIs{vv},sum(maskedVol(:)),contrasts{cc});

                if sum(maskedVol(:))>minVox
                    fprintf('%s %s %s: found %d voxels for %s\n',subinit,hemis{hh}, ROIs{vv},sum(maskedVol(:)),copenames{cc});
                    newVOI = thisVOI;
                    newVOI.vol = maskedVol;
                    savename = [VOIdir2 hemis{hh} '_' ROIs_4saving{vv} '_' copenames{cc} '.nii.gz'];
                    err = save_nifti(newVOI,savename);
                    if err~=0
                        error('error saving new VOI file')
                    end
                end
            end

        end

    end
    
end
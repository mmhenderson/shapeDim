%% Transform average brain (MNI average) from Talairach space to this 
% subject's anatomical scan space (and then into their functional space)

clear
close all

sub = 'S01';
sub_doreti = 'CG';
exp_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';

% this is the nifti defining MD network in TAL space
tal_parcels_orig = fullfile(exp_path, 'Processing/MDROI_XFM_TAL.nii.gz');

% this is the average tal brain
tal_brain_orig_mgz = '/mnt/neurocube/local/freesurfer6/average/mni305.cor.mgz';
tal_brain_orig_nii = fullfile(exp_path,'DataPreproc',sub,'mni305.cor.nii.gz');
if ~exist(tal_brain_orig_nii,'file')
    err = unix(['mri_convert ' tal_brain_orig_mgz ' ' tal_brain_orig_nii]);
    if err
        error('your mri_convert command failed!')
    end
end

% this is the MNI-average brain, after registering to subject's space
% ( we are about to make this)
tal_brain_xfm_nii = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_XFM.nii.gz');

% this is the MD-network parcellation, after registering to the subject's
% space
% (we are about to make this also)
tal_parcels_xfm = fullfile(exp_path,'DataPreproc',sub,'MDROI_XFM.nii.gz');

% write out file names for the reg files we will create now
manual_reg_file_dat = fullfile(exp_path,'DataPreproc',sub,'registration_manual_tal2anat.dat');
manual_reg_file_mat = fullfile(exp_path,'DataPreproc',sub,'registration_manual_tal2anat.mat');
auto_reg_file_dat = fullfile(exp_path,'DataPreproc',sub,'registration_auto_tal2anat.dat');
auto_reg_file_mat = fullfile(exp_path,'DataPreproc',sub,'registration_auto_tal2anat.mat');

% this is the subject's MC template from this experiment, need to use as a
% reference volume (made during main preproc script)
mctemplate = fullfile(exp_path,'DataPreproc',sub,'MCTemplate01.nii.gz');
% header info in this file - we also need this to make the transformations
% work correctly 
regheadercenter_file = fullfile(exp_path,'DataPreproc',sub,'regheadercenterMOD.mat');

%% manual registration to subject's anat
[s,w] = unix(['tkregister2 --s ' sub_doreti ' --mov ' tal_brain_orig_nii ' --regheader-center --reg ' manual_reg_file_dat]);
if s == 1 %a non-zero value returned by unix means the command has failed
    fprintf(['\n\nERROR: Failed to start manual registration. \n\n' 'ERROR MESSAGE FROM UNIX WRAPPER:' w '\n\n']);
    return
end

%% automatic fine-tuning
% (i am actually not using this because it seems to work poorly)
bbout_file = fullfile(exp_path,'DataPreproc',sub,'bbout_tal2anat.nii.gz');
fprintf('\n...Performing automatic registration using boundary based alignment, give it a minute...\n');
unix(['bbregister --s ' sub_doreti ' --mov ' tal_brain_orig_nii ' --reg ',...
    auto_reg_file_dat ' --init-reg ' manual_reg_file_dat ' --bold --o ' bbout_file]);
fprintf('Finished automatic registration!\n');

%% get ready to do the real transformation

% % convert the .dat file that we save out into a .mat file
% unix(['tkregister2 --s ', sub_doreti, ' --mov ', tal_brain_orig_nii ' --reg ',...
%     auto_reg_file_dat ' --fslregout ', auto_reg_file_mat ' --noedit']);

% convert the .dat file that we save out into a .mat file
unix(['tkregister2 --s ', sub_doreti, ' --mov ', tal_brain_orig_nii ' --reg ',...
    manual_reg_file_dat ' --fslregout ', manual_reg_file_mat ' --noedit']);
%%
% contatenate the header info and the registration
% name for the concat file
concat_file = fullfile(exp_path,'DataPreproc',sub,'concat_tal2anat.mat');
% unix(['convert_xfm -omat ' concat_file ' -concat ' regheadercenter_file ' ',...
%     auto_reg_file_mat]); 
unix(['convert_xfm -omat ' concat_file ' -concat ' regheadercenter_file ' ',...
    manual_reg_file_mat]); 

%% do the transformation of interest here

% this works roughly (as well as your registration worked)
err = unix(['flirt -in ', tal_brain_orig_nii, ' -ref ', mctemplate,...
    ' -out ', tal_brain_xfm_nii,' -init ', concat_file ' -applyxfm']);

err = unix(['flirt -in ', tal_parcels_orig, ' -ref ', mctemplate,...
    ' -out ', tal_parcels_xfm,' -init ', concat_file ' -applyxfm']);

%% fix the parcel labels so they don't have intermediate values

parcels_xfm_nii = load_nifti(tal_parcels_xfm);
vol = parcels_xfm_nii.vol;
labels = unique(vol(:));

mean_z = sum(sum(vol==labels(9), 2),1);
z2plot_mask = find(mean_z==max(mean_z));

figure;imagesc(vol(:,:,z2plot_mask));colorbar();title('z-slice')


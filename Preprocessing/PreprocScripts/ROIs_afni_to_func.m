function ROIs_afni_to_func()

% convert rois into functional space
% for the subjects where retinotopy was done in new afni pipeline
% need to run doreti_ROI2Vol (version in this folder) first

% note for this to work, you have to run getVisualROIs.m to make the
% Func2Anat_auto file. that script will crash for these subjects since rois
% are in the wrong format, but before it crashes it makes a file we need. 

cur_dir = pwd;

% ss = 8;
% SUBJ = 'DA';

ss = 9;
SUBJ = 'DF';


DORETI_path = '/mnt/neurocube/local/serenceslab/Doreti';
% NII_path = fullfile(DORETI_path,'RAW',SUBJ,'NII');
% ANAT_path = fullfile(DORETI_path,'ANAT',SUBJ);
% FUNC_path = fullfile(DORETI_path,'FUNC',SUBJ);
% addpath(genpath(fullfile(DORETI_path,'SCRIPTS_DOTS')));
AFNI_path = fullfile(DORETI_path,'FUNC',SUBJ,'PreProc','AFNI/'); % needs '/' to run smoothly!

roi_folder_in = sprintf('/usr/local/serenceslab/maggie/shapeDim/VOIs/from_afni/%s', SUBJ);
output_folder = sprintf('/usr/local/serenceslab/maggie/shapeDim/VOIs/S%02d/', ss);

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

%% get transformation file in correct format

preproc_folder = sprintf('/usr/local/serenceslab/maggie/shapeDim/DataPreproc/S%02d', ss);
xfm_file = fullfile(preproc_folder, 'Func2Anat_auto.dat');
xfm_file_mat = fullfile(preproc_folder, 'Func2Anat_auto.mat');
xfm_file_mat_inv = fullfile(preproc_folder, 'Func2Anat_auto_inverse.mat');
mctemplate_file = fullfile(preproc_folder, 'MCTemplateXFM01.nii.gz');

% convert to FSL format (.mat)
cmd = sprintf('tkregister2 --s %s --mov %s --reg %s --fslregout %s --noedit', ...
    SUBJ, mctemplate_file, xfm_file, xfm_file_mat);
fprintf(cmd)
fprintf('\n')
[s, ~] = unix(cmd);

% make the inverse transfrmation matrix (anat to func)
cmd = sprintf('convert_xfm -omat %s -inverse %s', ...
    xfm_file_mat_inv, xfm_file_mat);
fprintf(cmd)
fprintf('\n')
[s, ~] = unix(cmd);

%% apply transform from anat to func space

% first I am doing this to my actual anat volume, just to check the 
% alignment is correct. It's hard to tell with the ROI volumes. But the 
% ROIs I have are in same space as high-res anatomy, so if this
% registration looks right, then the ROIs should be correct.

anatfile = fullfile(DORETI_path, 'ANAT',SUBJ, 'mri', 'orig.nii.gz');
checkfile = fullfile(output_folder, 'anat_reg2func_check.nii.gz');
% ^ check this file against MCTemplate in fsleyes

cmd = sprintf('flirt -in %s -ref %s -out %s -init %s -applyxfm', ...
    anatfile, mctemplate_file, checkfile, xfm_file_mat_inv);

fprintf(cmd)
fprintf('\n')
[s, ~] = unix(cmd);

%% now transform ROIs

my_rois = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','V3AB','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
% my_rois = {'hV4','LO1','LO2'};

hemis = {'lh','rh'};
for hemi = 1:2
    for rr = 1:numel(my_rois)
        roi_name = [hemis{hemi} '_' my_rois{rr}]; % base name for roi labels
        
        roi_filename_orig = fullfile(roi_folder_in, roi_name);
        roi_filename_xfm = fullfile(output_folder, [roi_name, '_xfm.nii.gz']);
        roi_filename_xfm_thresh = fullfile(output_folder,  [roi_name, '_xfm_thresh.nii.gz']);
        roi_filename_binary = fullfile(output_folder,  [roi_name, '.nii.gz']);
        
        % first transform the file
        cmd = sprintf('flirt -in %s -ref %s -out %s -init %s -applyxfm', ...
            roi_filename_orig, mctemplate_file, roi_filename_xfm, xfm_file_mat_inv);

        fprintf(cmd)
        fprintf('\n')
        [s, ~] = unix(cmd);
        
        % then apply threshold to make it binary again
        % setting threshold of 0.3 here, this is what was in old code for
        % call to mris_label2vol
        thr = 0.3;
        cmd = sprintf('fslmaths %s -thr %.1f %s', ...
            roi_filename_xfm, thr, roi_filename_xfm_thresh);

        fprintf(cmd)
        fprintf('\n')
        [s, ~] = unix(cmd);
        
        % convert to binary 0/1
        cmd = sprintf('fslmaths %s -bin %s', ...
            roi_filename_xfm_thresh, roi_filename_binary);

        fprintf(cmd)
        fprintf('\n')
        [s, ~] = unix(cmd);
    end
end

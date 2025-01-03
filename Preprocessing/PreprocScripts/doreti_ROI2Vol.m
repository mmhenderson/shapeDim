function doreti_ROI2Vol(SUBJ)

% MMH copied this out of the doreti folder so that i can edit and use as needed for shapedim
% To be used for subjects DA, DF, that had retinotopy drawn in AFNI format.
% we need to convert the roi files into a format that will work with fsl/freesurfer
% this is step 1, then run ROIs_afni_to_func next

cur_dir = pwd;

DORETI_path = '/mnt/neurocube/local/serenceslab/Doreti';
% NII_path = fullfile(DORETI_path,'RAW',SUBJ,'NII');
% ANAT_path = fullfile(DORETI_path,'ANAT',SUBJ);
% FUNC_path = fullfile(DORETI_path,'FUNC',SUBJ);
% addpath(genpath(fullfile(DORETI_path,'SCRIPTS_DOTS')));
AFNI_path = fullfile(DORETI_path,'FUNC',SUBJ,'PreProc','AFNI/'); % needs '/' to run smoothly!

output_folder = sprintf('/usr/local/serenceslab/maggie/shapeDim/VOIs/from_afni/%s', SUBJ);
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

%% POST-ROI-DRAWING
%cd(AFNI_path)
my_rois = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','V3AB','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
% my_rois = {'hV4','LO1','LO2'};
% my_rois = {'V1d','V2d','V3d'};
% my_rois = {'MT'};
hemis = {'lh','rh'};
for hemi = 1:2
    for rr = 1:numel(my_rois)
        roi_name = [hemis{hemi} '_' my_rois{rr}]; % base name for roi labels
        
        roi_filename = fullfile(AFNI_path, roi_name);
        dset_filename = fullfile(output_folder, roi_name);
        nii_filename = fullfile(output_folder, roi_name);
        
%         grid_parent_filename = fullfile(AFNI_path, 'REG_MC_DET_ZNORM_VB_avg.nii.gz');
        grid_parent_filename = fullfile(DORETI_path, 'ANAT', SUBJ, 'mri','orig.nii.gz');
        % ^ so the ROI volumes should end up in same space as high res
        % anatomy
        
        % convert .roi files to .dset files
        unix(['ROI2dataset -prefix ' dset_filename ' \'...
            '-of 1D \'...
            '-input ' roi_filename '.1D.roi']);
        
        % convert .dset files (surface) to .nii.gz files (volume)
        [s,~] = unix(['3dSurf2Vol -spec $SUBJECTS_DIR/' SUBJ '/SUMA/' SUBJ '_' hemis{hemi} '.spec \'...
            '-surf_A smoothwm \'...
            '-surf_B pial \'...
            '-sv $SUBJECTS_DIR/' SUBJ '/SUMA/' SUBJ '_SurfVol.nii \'...
            '-grid_parent ' grid_parent_filename ' \'...
            '-map_func max \'...
            '-f_steps 10 \'...
            '-f_p1_fr      -0.1  \'...
            '-f_pn_fr      -0.1  \'...
            '-f_index voxels  \'...
            '-sdata_1D ' dset_filename '.1D.dset \'...
            '-prefix ' nii_filename '.nii.gz']);

    end
end

end

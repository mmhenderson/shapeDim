%% Transform average brain (MNI average) from Talairach space to this 
% subject's anatomical scan space (and then into their functional space)

clear
close all

root = pwd;
filesepinds = find(root==filesep);
exp_path = root(1:filesepinds(end));
subinit = 'CG';
subnum = 'S01';
subretfolder = fullfile(root(1:filesepinds(end-3)), 'Doreti','ANAT',subinit);

% this is the nifti defining MD network in TAL space
tal_parcels_orig = fullfile(exp_path, 'Processing/MDROI_XFM_TAL.nii.gz');

% this is the average tal brain
tal_brain_orig_nii = fullfile(exp_path,'DataPreproc',subnum,'mni305.cor.nii.gz');

% this is the MC template (final space we want to end up in
mctemplate = fullfile(exp_path,'DataPreproc',subnum,'MCTemplateXFM01.nii.gz');





%% load the average brain 

avg_brain_nii = load_nifti(tal_brain_orig_nii);

% get out some important fields from it
sform = avg_brain_nii.sform;
vol = avg_brain_nii.vol;
dims = size(vol);

%% load the parcellation of MD regions 

md_parcels_nii = load_nifti(tal_parcels_orig);

rois_vol = md_parcels_nii.vol;

labels = unique(rois_vol(:));
%% create a mask of the original brain 
% this is only for visualization/making sure the transformation worked,
% since it's hard to tell if the MD ROIs are lined up w brain
% first get all the coordinates of the points above some arbitrary threshold 
mask = double(vol>100);
[i,j,k] = ind2sub(size(mask),find(mask==1));
% here i'm converting the points to mm, which we will need down below when
% we do the conversion into subject anatoical space. 
mask_mm = coord2mm([i,j,k],sform);  
% converting back to coordinates (matrix indices) again - these should be
% same as i,j,k above.
mask_coords = mm2coord(mask_mm,sform);
% then back to linear indices 
mask_inds = sub2ind(dims, mask_coords(:,1), mask_coords(:,2), mask_coords(:,3));

mask_new = zeros(size(mask));
mask_new(mask_inds) =1;

% make and save the new nifti
mask_nii = avg_brain_nii;
mask_nii.vol = mask_new;
tal_brain_mask_nii = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK.nii.gz');
save_nifti(mask_nii,tal_brain_mask_nii);

%% load the talairach-anat transformation matrix
% load this file, created by recon-all
tal_reg_file = fullfile(subretfolder, 'mri','transforms', 'talairach.xfm');

% read it to get my tranformation matrix
fid = fopen(tal_reg_file,'r');
tal_xfm_mat = [];
line = fgetl(fid);
count = 0;
while ischar(line) 
    count = count+1;
    line = fgetl(fid);
    if count>=8
        spaces = [0,find(line==' ')];
        vals = [];
        for sp =1:length(spaces)-1
            vals = [vals, str2double(line(spaces(sp)+1:spaces(sp+1)-1))];
        end
        tal_xfm_mat = [tal_xfm_mat; vals];
    end        
end
fclose(fid);

tal_rot_scale = tal_xfm_mat(1:3,1:3);
tal_translate = tal_xfm_mat(1:3,4)';

%% convert the mask into anatomical space.

% apply the transformation matrix
mask_mm_xfm = mask_mm*pinv(tal_rot_scale) - tal_translate;

% get from mm back to matrix indices
mask_coords_xfm = round(mm2coord(mask_mm_xfm,sform));

% make sure no indices are out of range
out_of_bounds = any(mask_coords_xfm>dims,2) |...
    any(mask_coords_xfm<1);
mask_coords_xfm(out_of_bounds) = [];

% back to linear indices
mask_inds_xfm = sub2ind(dims, mask_coords_xfm(:,1), mask_coords_xfm(:,2), mask_coords_xfm(:,3));
mask_xfm = zeros(dims);
mask_xfm(mask_inds_xfm) = 1;

% make and save the new nifti
mask_xfm_nii = avg_brain_nii;
mask_xfm_nii.vol = mask_xfm;
tal_brain_mask_xfm = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK_XFMANAT.nii.gz');
save_nifti(mask_xfm_nii,tal_brain_mask_xfm);

%% convert the MD ROIs into anatomical space.

roi_vol_xfm = zeros(size(rois_vol));

for ll=1:length(labels)

    % get all coordinates with this label
    [i,j,k] = ind2sub(dims, find(rois_vol==labels(ll)));
    roi_mm = coord2mm([i,j,k],sform);
    
    % apply the transformation
    roi_mm_xfm = roi_mm*pinv(rot_scale) - translate;
    
    % get from mm back to matrix indices
    roi_coords_xfm = round(mm2coord(roi_mm_xfm,sform));
    
    % make sure no indices are out of range
    out_of_bounds = any(roi_coords_xfm>dims,2) |...
        any(roi_coords_xfm<1,2);
    roi_coords_xfm(out_of_bounds,:) = [];

    roi_inds_xfm = sub2ind(dims, roi_coords_xfm(:,1), roi_coords_xfm(:,2), roi_coords_xfm(:,3));

    roi_vol_xfm(roi_inds_xfm) = labels(ll);
    
end

% make and save the nifti
md_xfm_nii = avg_brain_nii;
md_xfm_nii.vol = roi_vol_xfm;
md_xfm_file = fullfile(exp_path,'DataPreproc', subnum, 'MDROI_XFMANAT.nii.gz');
save_nifti(md_xfm_nii,md_xfm_file);

%% load the mc template

mctemplate_nii = load_nifti(mctemplate);
mctemplate_vol = mctemplate_nii.vol;
mctemplate_dims = size(mctemplate_vol);
mctemplate_sform = mctemplate_nii.sform;

%% load the anat to funct registration
func2anat_reg_dat = fullfile(exp_path,'DataPreproc',subnum,'Func2Anat_auto.dat');
func2anat_reg_mat = fullfile(exp_path,'DataPreproc',subnum,'Func2Anat_auto.mat');

if ~exist(func2anat_reg_mat,'file')
    % convert the .dat file that we save out into a .mat file
    unix(['tkregister2 --s ', subinit, ' --mov ', mctemplate ' --reg ',...
       func2anat_reg_dat ' --fslregout ', func2anat_reg_mat ' --noedit']);
end

% read it to get my tranformation matrix
fid = fopen(func2anat_reg_mat,'r');
xfm_mat = [];
line = 'x';
count = 0;
while ischar(line) && size(xfm_mat,1)<3
    count = count+1;
    line = fgetl(fid);
    spaces = [0,find(line==' ')];
    vals = [];
    for sp =1:length(spaces)-1
        newval =str2double(line(spaces(sp)+1:spaces(sp+1)-1));
        if ~isnan(newval)
            vals = [vals, newval];
        end
    end
    xfm_mat = [xfm_mat; vals];
end
fclose(fid);
xfm_mat = [xfm_mat; [0,0,0,1]];
rot_scale = xfm_mat(1:3,1:3);
translate = xfm_mat(1:3,4)';

%% load the mask from anat space

tal_brain_mask_xfm = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK_XFMANAT.nii.gz');
mask_anat = load_nifti(tal_brain_mask_xfm);
anat_mask_vol = mask_anat.vol;
anat_dims = size(anat_mask_vol);
anat_sform = mask_anat.sform;

[i,j,k] = ind2sub(anat_dims, find(anat_mask_vol==1));
mask_coords = [i,j,k];
mask_mm = coord2mm([i,j,k], anat_sform);

%% convert the mask into functional space

coord_orig = [mask_coords(1,:)+1,1];

full_xfm_mat = pinv(anat_sform)*pinv(xfm_mat)*pinv(mctemplate_sform);

coord_trans = round(coord_orig*full_xfm_mat)

%%
% mask_coords_xfm = round(mask_mm_corner*pinv(rot_scale) - translate) + corner_mm;

% get from mm back to matrix indices (in the new mctemplate matrix)
% mask_coords_xfm = round(mm2coord(mask_mm_xfm,mctemplate_sform));

mask_coords_xfm(1:10,:)

% make sure no indices are out of range
out_of_bounds = any(mask_coords_xfm>mctemplate_dims,2) |...
    any(mask_coords_xfm<1,2);
mask_coords_xfm(out_of_bounds,:) = [];

% back to linear indices
mask_inds_xfm = sub2ind(mctemplate_dims, mask_coords_xfm(:,1), mask_coords_xfm(:,2), mask_coords_xfm(:,3));
mask_xfm = zeros(mctemplate_dims);
mask_xfm(mask_inds_xfm) = 1;

% make and save the new nifti
mask_xfm_nii = mctemplate_nii;
mask_xfm_nii.vol = mask_xfm;
tal_brain_mask_xfm = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK_XFMFUNC.nii.gz');
save_nifti(mask_xfm_nii,tal_brain_mask_xfm);

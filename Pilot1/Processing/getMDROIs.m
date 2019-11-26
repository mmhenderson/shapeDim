%% Create MD network ROI definitions for individual subjects
% first project from Talairach space to anatomical, then from anatomical
% space into functional space. finally, interesect the definitions in
% functional space with the activated regions from MD localizer task.

% MMH 11/26/19
%%

clear
close all

root = pwd;
filesepinds = find(root==filesep);
exp_path = root(1:filesepinds(end));
subinit = 'CG';
subnum = 'S01';
subretfolder = fullfile(root(1:filesepinds(end-3)), 'Doreti','ANAT',subinit);

% this is the nifti defining MD network in TAL space (make this using
% getMD_avg.m)
md_parcels_orig_nii = fullfile(exp_path, 'Processing/MDROI_XFM_TAL.nii.gz');

% this is the average MNI brain (there is a copy of this at 
% '/mnt/neurocube/local/freesurfer6/average/mni305.cor.mgz' 
% and should also be one in this folder.
avg_brain_orig_nii = fullfile(exp_path,'DataPreproc',subnum,'mni305.cor.nii.gz');

% make the roi directory
md_roi_dir = fullfile(exp_path,'VOIs_MD',subnum);
if ~isfolder(md_roi_dir)
    mkdir(md_roi_dir);
end

main_roi_dir = fullfile(exp_path,'VOIs',subnum);

% need this file to get from anat into functional space
regheadercenter_file = fullfile(exp_path,'DataPreproc',subnum,'regheadercenterMOD.mat');

% need this file for a template when converting into functional space
mctemplate_file = fullfile(exp_path,'DataPreproc',subnum, 'MCTemplateXFM01.nii.gz');

% make just these VOIs (not every single parcel)
VOIs = {'IFS', 'AI-FO', 'iPCS', 'sPCS','poCS'};
nVOIs = length(VOIs);

% list which parcels in the MD file correspond to which VOIs in the
% numbered list (numbered roughly anterior:posterior, column 1 is left, 2 is right)
parcellist = {  [23,24],[21,22];
                [9],[8];
                [11],[10];
                [14],[13];
                [5],[7]};

hemis = {'lh','rh'};
nHemis = length(hemis);

% minimum number of voxels to count it as a ROI
minvox = 10;
%% load the average brain 
avg_brain = load_nifti(avg_brain_orig_nii);
% sform tells how to get from voxel coordinates in the matrix, into mm
% coordinates in standard space.
orig_sform = avg_brain.sform;
avg_brain_vol = avg_brain.vol;
orig_dims = size(avg_brain_vol);

%% load the parcellation of MD regions 

md_parcels = load_nifti(md_parcels_orig_nii);
rois_vol = md_parcels.vol;
labels = unique(rois_vol(:));
% this file has same dims and sform as the average brain.

%% load subject's anatomy file (use as a template for first transformation)
template_nii = fullfile(subretfolder,'mri','orig.nii.gz');
template = load_nifti(template_nii);

% get out some important fields from it
template_sform = template.sform;
template_vol = template.vol;
template_dims = size(template_vol);

%% create a mask of the original MNI average brain

% this is only for visualization/making sure the transformation worked,
% since it's hard to tell if the MD ROIs are lined up w brain.
% first get all the coordinates of the points above some arbitrary threshold 
avg_brain_mask = double(avg_brain_vol>100);
[i,j,k] = ind2sub(size(avg_brain_mask),find(avg_brain_mask==1));

% convert these voxel coordinates into mm, using the sform matrix
mask_mm = coord2mm([i,j,k],orig_sform); 

% now convert from mm to voxel coords in the NEW matrix, using the sform
% matrix from the template file. This way it will be in the same coordinate
% system as the subject's anatomy file.
mask_coords = round(mm2coord(mask_mm,template_sform));

% make sure no indices are out of range
out_of_bounds = any(mask_coords>template_dims,2) |...
    any(mask_coords<1,2);
mask_coords(out_of_bounds,:) = [];

% then back to linear indices 
mask_inds = sub2ind(template_dims, mask_coords(:,1), mask_coords(:,2), mask_coords(:,3));

mask = zeros(template_dims);
mask(mask_inds) =1;

% make and save the new nifti (this has NOT been transformed, but has been
% switched to a new coordinate system). 
mask_new = template;
mask_new.vol = mask;
tal_brain_mask_nii = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK.nii.gz');
save_nifti(mask_new,tal_brain_mask_nii);
fprintf('saved mask file to %s\n',tal_brain_mask_nii);
%% load the talairach-anat transformation matrix
% NOTE that this is not a rigid-body transform, it is a 12 parameter affine
% transformation so it includes scaling etc.

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

% these map from anat space to tal space, so we will have to apply them in reverse
% in order to get from tal space to anat
rot_scale = tal_xfm_mat(1:3,1:3);
translate = tal_xfm_mat(1:3,4)';

%% convert the mask into anatomical space.
% this can be used to double check and make sure the alignment works, but
% is not strictly necessary

% apply the transformation matrix
mask_mm_xfm = mask_mm*pinv(rot_scale) - translate;

% get from mm back to matrix indices
mask_coords_xfm = round(mm2coord(mask_mm_xfm,template_sform));

% make sure no indices are out of range
out_of_bounds = any(mask_coords_xfm>template_dims,2) |...
    any(mask_coords_xfm<1,2);
mask_coords_xfm(out_of_bounds,:) = [];

% back to linear indices
mask_inds_xfm = sub2ind(template_dims, mask_coords_xfm(:,1), mask_coords_xfm(:,2), mask_coords_xfm(:,3));
mask_xfm = zeros(template_dims);
mask_xfm(mask_inds_xfm) = 1;

% make and save the new nifti
mask_xfm_struct = template;
mask_xfm_struct.vol = mask_xfm;
tal_brain_mask_xfm_nii = fullfile(exp_path,'DataPreproc', subnum, 'MNI_BRAIN_MASK_REG2ANAT.nii.gz');
save_nifti(mask_xfm_struct,tal_brain_mask_xfm_nii);
fprintf('saved transformed (anat) mask file to %s\n',tal_brain_mask_xfm_nii);
%% convert the MD ROIs into anatomical space.

% making a big mask that has all the ROIs
rois_all_xfm = zeros(size(rois_vol));

for ll=1:length(labels)

    % get all coordinates with this label
    [i,j,k] = ind2sub(orig_dims, find(rois_vol==labels(ll)));
    roi_mm = coord2mm([i,j,k],orig_sform);
    
    % apply the transformation
    roi_mm_xfm = roi_mm*pinv(rot_scale) - translate;
    
    % get from mm back to matrix indices (now indexing into template)
    roi_coords_xfm = round(mm2coord(roi_mm_xfm,template_sform));
    
    % make sure no indices are out of range
    out_of_bounds = any(roi_coords_xfm>template_dims,2) |...
        any(roi_coords_xfm<1,2);
    roi_coords_xfm(out_of_bounds,:) = [];

    % back to linear indices
    roi_inds_xfm = sub2ind(template_dims, roi_coords_xfm(:,1), roi_coords_xfm(:,2), roi_coords_xfm(:,3));

    % put the label into big matrix
    rois_all_xfm(roi_inds_xfm) = labels(ll);
    
    % also making and saving a nifti with the single ROI label (to be used
    % for making ROIs later on)
    roi_single_xfm = zeros(size(rois_vol));
    roi_single_xfm(roi_inds_xfm) = 1;
    
    roi_xfm = template;
    roi_xfm.vol = roi_single_xfm;
    roi_xfm_nii = fullfile(md_roi_dir,sprintf('MDROI_label%d_REG2ANAT.nii.gz',labels(ll)));
    save_nifti(roi_xfm, roi_xfm_nii);
    fprintf('saved MD ROI (anat) label %d file to %s\n',ll,roi_xfm_nii);
end

% make and save the (full) nifti
md_xfm = template;
md_xfm.vol = rois_all_xfm;
md_xfm_nii = fullfile(exp_path,'DataPreproc', subnum, 'MDROI_REG2ANAT.nii.gz');
save_nifti(md_xfm,md_xfm_nii);
fprintf('saved all MD ROIs (anat) to %s\n',md_xfm_nii);

%% convert all ROIs from anatomical to functional space

for ll = 1:length(labels)

    % this is the thing we just made
    file2align = fullfile(md_roi_dir,sprintf('MDROI_label%d_REG2ANAT.nii.gz',labels(ll)));
    
    % this is the transformed ROI
    file2save = fullfile(md_roi_dir,sprintf('MDROI_label%d_REG2FUNC.nii.gz',labels(ll)));
   
    % make the transformation w flirt
    err =  unix(['flirt -in ', file2align ' -ref ', mctemplate_file,...
        ' -applyxfm -init ',regheadercenter_file ' -out ', file2save]);
 
    fprintf('saved MD ROI (functional) label %d file to %s\n',ll,file2save);
    
	% apply a threshold to re-binarize the mask 
    % (the flirt output is interpolated so it has continuous values 0-1)
    new_nii = load_nifti(file2save);
    new_nii.vol = double(new_nii.vol>0.5);  % 0.5 will keep the mask around the same size
    save_nifti(new_nii,file2save);
end

%% load MD Localizer t-map
currdir = pwd;
% cope 3 is hard>easy
stats_path = fullfile(exp_path,'AnalyzeMDLocalizer',subnum,'feats','AllSessionsFE.gfeat','cope3.feat','stats');
cd(stats_path)
% load my t-stat map for all voxels
nifti = load_nifti('tstat1.nii.gz');
vmpVals = nifti.vol; %Store t-values for each voxel index that is in my visual areas in 'vmpVals'

clear nifti

% find my threshold 
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

% define my final mask
vals2use = vmpVals>my_t;

cd(currdir)
%% now intersect this with each ROI in funct space

for vv=1:nVOIs
        
    
        for hh=1:nHemis
            parcels = parcellist{vv,hh};
    
            
            % mask out this ROI
            for pp=1:length(parcels)

                label_full_nii = fullfile(md_roi_dir,sprintf('MDROI_label%d_REG2FUNC.nii.gz',parcels(pp)));
               
                label_full = load_nifti(label_full_nii);
                                 
                if pp==1
                    fullVol = zeros(size(label_full.vol));
                end
                
                fullVol = fullVol | label_full.vol;
                
            end

            maskedVol = fullVol & vals2use;

            if sum(maskedVol(:))>minvox
                fprintf('%s %s %s: found %d voxels for hard>easy\n',subinit,hemis{hh}, VOIs{vv},sum(maskedVol(:)));
                newVOI = label_full;
                newVOI.vol = maskedVol;
                savename = fullfile(main_roi_dir,sprintf('%s_%s_mdloc.nii.gz',hemis{hh},VOIs{vv}));

                err = save_nifti(newVOI,savename);
                if err~=0
                    error('error saving new VOI file')
                end
            end
          
            
        end
        
end

%% Transform average brain (MNI average) from Talairach space to this 
% subject's anatomical scan space (and then into their functional space)

clear
close all

sub = 'S01';
sub_doreti = 'CG';
exp_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';

subretfolder = fullfile('/mnt/neurocube/local/serenceslab/Doreti/','ANAT',sub_doreti);

% this is the nifti defining MD network in TAL space
tal_parcels_orig = fullfile(exp_path, 'Processing/MDROI_XFM_TAL.nii.gz');

% this is the average tal brain
tal_brain_orig_nii = fullfile(exp_path,'DataPreproc',sub,'mni305.cor.nii.gz');

% this is the subject's MC template from this experiment, need to use as a
% reference volume (made during main preproc script)
mctemplate = fullfile(exp_path,'DataPreproc',sub,'MCTemplate01.nii.gz');

%% this is the concat file for the transformation, already made it

% without the header
manual_reg_file_mat = fullfile(exp_path,'DataPreproc',sub,'registration_manual_tal2anat.mat');
% contatenate the header info and the registration
% name for the concat file
concat_file = fullfile(exp_path,'DataPreproc',sub,'concat_tal2anat.mat');

%% load the average brain 

avg_brain_nii = load_nifti(tal_brain_orig_nii);

sform = avg_brain_nii.sform;
vol = avg_brain_nii.vol;
dims = size(vol);

% these are the actual mm coords (in original brain) that correspond to
% each array in the index.
vox_mm_list = [(1:dims(1))',(1:dims(2))',(1:dims(3))',ones(dims(1),1)]*sform';
vox_mm_list = vox_mm_list(:,1:3);
x_mm = vox_mm_list(:,1);
y_mm = vox_mm_list(:,2);
z_mm = vox_mm_list(:,3);

%% create a mask of the original brain
mask = double(vol>100);
avg_brain_nii.vol = mask;
tal_brain_mask_nii = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK.nii.gz');
save_nifti(avg_brain_nii,tal_brain_mask_nii);

%% make another mask, to check understanding of this coordinate system
% (this one is identical to the above mask)
[i,j,k] = ind2sub(size(mask),find(mask==1));
mask_mm = coord2mm([i,j,k],sform);
mask_coords = mm2coord(mask_mm,sform);
mask_inds = sub2ind(dims, mask_coords(:,1), mask_coords(:,2), mask_coords(:,3));

mask_new = zeros(size(mask));
mask_new(mask_inds) =1;

avg_brain_nii.vol = mask_new;
tal_brain_mask_nii2 = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK2.nii.gz');
save_nifti(avg_brain_nii,tal_brain_mask_nii2);

%% load the transformation matrix

% read it to get my tranformation matrix
fid = fopen(manual_reg_file_mat,'r');
xfm_mat_manual = [];
line = 'x';
count = 0;
while ischar(line) 
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
    xfm_mat_manual = [xfm_mat_manual; vals];
end
fclose(fid);

% read the concat file too
fid = fopen(concat_file,'r');
xfm_mat_concat = [];
line = 'x';
count = 0;
while ischar(line) 
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
    xfm_mat_concat = [xfm_mat_concat; vals];
end
fclose(fid);

regheadercentermod = [232, 220, 210];
xfm_mat = xfm_mat_concat;
xfm_mat(1:3,4) = xfm_mat(1:3,4) - regheadercentermod';

xfm_mat = xfm_mat_manual;

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

%% convert the mask with flirt
tal_brain_mask_xfm1 = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK_XFM1.nii.gz');
err = unix(['flirt -in ', tal_brain_mask_nii, ' -ref ', mctemplate,...
    ' -out ', tal_brain_mask_xfm1,' -init ', concat_file ' -applyxfm']);

%% load mctemplate to get its sform matrix
mctemplate_nii = load_nifti(mctemplate);
new_sform = mctemplate_nii.sform;
new_dims = size(mctemplate_nii.vol);

%% put the mask into the new coord system, without transforming 
% this is anoher check to make sure we understand the coordinate systems
[i,j,k] = ind2sub(size(mask),find(mask==1));
mask_mm = coord2mm([i,j,k],sform);
mask_coords_new = mm2coord(mask_mm,new_sform);
for dd=1:3
    mask_coords_new(:,dd) = max(min(mask_coords_new(:,dd),new_dims(dd)),1);
end
mask_coords_new = round(mask_coords_new);
mask_inds = sub2ind(new_dims, mask_coords_new(:,1), mask_coords_new(:,2), mask_coords_new(:,3));

mask_new = zeros(new_dims);
mask_new(mask_inds) =1;

mctemplate_nii.vol = mask_new;
tal_brain_mask_nii3 = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK3.nii.gz');
save_nifti(mctemplate_nii,tal_brain_mask_nii3);

%% now convert the mask manually (to funct space)
rot_scale = xfm_mat(1:3,1:3);
translate = xfm_mat(1:3,4)';

mask_mm_xfm = mask_mm*rot_scale + translate;

new_sform_xfm = new_sform*xfm_mat;

mask_coords_xfm = round(mm2coord(mask_mm_xfm,new_sform_xfm));
for dd=1:3
    mask_coords_xfm(:,dd) = max(min(mask_coords_xfm(:,dd),new_dims(dd)),1);
end

mask_inds_xfm = sub2ind(new_dims, mask_coords_xfm(:,1), mask_coords_xfm(:,2), mask_coords_xfm(:,3));
mask_xfm = zeros(new_dims);
mask_xfm(mask_inds_xfm) = 1;

mctemplate_nii.vol = mask_xfm;
tal_brain_mask_xfm2 = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK_XFM2.nii.gz');
save_nifti(mctemplate_nii,tal_brain_mask_xfm2);

%% now convert the mask manually  (to anat space)
rot_scale = tal_xfm_mat(1:3,1:3);
translate = tal_xfm_mat(1:3,4)';

mask_mm_xfm = mask_mm*pinv(rot_scale) - translate;

mask_coords_xfm = round(mm2coord(mask_mm_xfm,sform));
for dd=1:3
    mask_coords_xfm(:,dd) = max(min(mask_coords_xfm(:,dd),dims(dd)),1);
end

mask_inds_xfm = sub2ind(dims, mask_coords_xfm(:,1), mask_coords_xfm(:,2), mask_coords_xfm(:,3));
mask_xfm = zeros(dims);
mask_xfm(mask_inds_xfm) = 1;

avg_brain_nii.vol = mask_xfm;
tal_brain_mask_xfm3 = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK_XFM3.nii.gz');
save_nifti(avg_brain_nii,tal_brain_mask_xfm3);


%%

% %% put a thing in it to test coordinates
% center = abs(sform(1:3,4))';
% 
% mm2edit = [-40,-50,-60];
% ind2edit = mm2coord(mm2edit, sform);
% 
% mask_edited = mask;
% mask_edited(ind2edit(1),ind2edit(2),ind2edit(3)) = 12;
% 
% avg_brain_nii.vol = mask_edited;
% nii_new = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK_EDITED.nii.gz');
% save_nifti(avg_brain_nii,nii_new);
% 

% %% this is after applying sform to flip the coords around
% % taking the coords and multiplying by transpose of sform
% % x-coord is negative of first original column (left-right)
% % y-coord is positive of third original column (post-ant)
% % z-coord is negative of second original column (inf-superior)
% close all
% mask_adj = flip(mask,1);
% mask_adj = flip(mask_adj,2);
% mask_adj = permute(mask_adj,[1,3,2]);
% figure;imagesc(mask_adj(:,:,100)');xlabel('X');ylabel('Y');axis square % (remember imagesc plots rows on the y-axis)
% figure;imagesc(squeeze(mask_adj(:,100,:))');xlabel('X');ylabel('Z');axis square
% figure;imagesc(squeeze(mask_adj(100,:,:))');xlabel('Y');ylabel('Z');axis square
%% transform the mask w coords

[i,j,k] = ind2sub(size(mask),find(mask==1));
mask_mm = coord2mm([i,j,k],sform);

mask_xfm = zeros(size(mask));
% put these mm coords into the brain now
for mm=1:length(mask_mm_coords)

    % find all voxels that are close to this point (don't leave any
    % holes!)         
    xc = find(abs(x_mm-mask_mm_coords(mm,1))<voxmm);
    yc = find(abs(y_mm-mask_mm_coords(mm,3))<voxmm);
    zc = find(abs(z_mm-mask_mm_coords(mm,2))<voxmm);

    % label all of those voxels with the same label they had before
    [xc,yc,zc] = meshgrid(xc,yc,zc);
    xc = xc(:);yc = yc(:);zc = zc(:);
    for ii=1:numel(xc) 
        mask(xc(ii),yc(ii),zc(ii)) = 1;
    end

end

%% save the new nifti

new_nii = avg_brain_nii;
new_nii.vol = mask;
tal_brain_mask_nii = fullfile(exp_path,'DataPreproc', sub, 'MNI_BRAIN_MASK.nii.gz');
save_nifti(new_nii,tal_brain_mask_nii);

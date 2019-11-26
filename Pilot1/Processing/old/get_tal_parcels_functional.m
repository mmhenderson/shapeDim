% Makes a nii file that defines MD network in taleirach space.

% load the anatomical definitions of a parcellation of the MD network 
% from Fedorenko, Duncan , & Kanwisher, 2013
% downloaded from imaging.mrc-cbu.cam.ac.uk/imaging/MDsystem. 
% This script will convert that nifti file into a new nifti in
% the same coordinate system of the MNI-average brain from freesurfer. 
% (all happens in TAL space).
% Next step is to project this file into subject (functional) space using
% alignment between MNI-average brain and subject anatomical.

% MMH 11/22/19
%% 
clear

root = pwd;
filesepinds = find(root==filesep);
exp_path = root(1:filesepinds(end));
subinit = 'CG';
subnum = 'S01';
subretfolder = fullfile(root(1:filesepinds(end-3)), 'Doreti','ANAT',subinit);

%% define the transformation file
% load this file, created by recon-all
tal_reg_file = fullfile(subretfolder, 'mri','transforms', 'talairach.xfm');

% read it to get my tranformation matrix
fid = fopen(tal_reg_file,'r');
xfm_mat = [];
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
        xfm_mat = [xfm_mat; vals];
    end        
end
fclose(fid);

%% load the nifti with the MD ROIs masked out
% need to have a copy of this file in our current folder, this is the thing
% we downloaded from Fedorenko paper
parcels_nii = load_nifti(fullfile(exp_path, 'Processing', 'MDROI.nii.gz'));

% this is the name of the new file we will make here
parcels_xfm_nii = fullfile(exp_path, 'Processing', 'MDROI_XFM_TAL.nii.gz');

% this is the name of the second file we will make, with areas defined in
% subject anatomical space.
parcels_anat_nii = fullfile(exp_path,'DataPreproc',subnum,'MDROI_XFM_ANAT.nii.gz');


% range of the image in mm
% got these by loading in fslview, can also read them out from the header
% probably.
% x_range = [-78, 78];
x_range = [-79,79];
% y_range = [-112, 76];
y_range = [-113 , 77];
% z_range = [-50, 86];
z_range=  [-51, 87];

% size of the voxels (also in header)
voxmm = 2;

% the coordinate of each pt is the center of this voxel in mm
xlist = (x_range(1)+1:voxmm:x_range(2))';
ylist = (y_range(1)+1:voxmm:y_range(2))';
zlist = (z_range(1)+1:voxmm:z_range(2))';

% different numbers for each region in the MD network
labels = unique(parcels_nii.vol(:));

% check dims
parcel_dims = size(parcels_nii.vol);
if any(parcel_dims~=[length(xlist),length(ylist),length(zlist)])
    error('nifti image is wrong size')
end

%% load the MNI average brain
% this file is same for all subjects, it's a reference brain
tal_brain_orig_mgz = '/mnt/neurocube/local/freesurfer6/average/mni305.cor.mgz';
% make sure we have a copy of it in this folder, in nii format
tal_brain_orig_nii = fullfile(exp_path,'Processing','mni305.cor.nii.gz');
if ~exist(tal_brain_orig_nii,'file')
    err = unix(['mri_convert ' tal_brain_orig_mgz ' ' tal_brain_orig_nii]);
    if err
        error('your mri_convert command failed!')
    end
end
avg_brain = load_nifti(tal_brain_orig_nii);
avg_brain_dims = size(avg_brain.vol);

% list the mm coordinates of every voxel in this array
% sform tells us a mapping from i,j,k, in the volume to x,y,z in mm
avg_brain_sform = avg_brain.sform;
[i,j,k] = meshgrid(1:avg_brain_dims(1),1:avg_brain_dims(2),1:avg_brain_dims(3));
avg_brain_mm_list = [i(:),j(:),k(:), ones(numel(i(:)),1)]*avg_brain_sform';
avg_brain_mm_list = avg_brain_mm_list(:,1:3);

% taking unique in each dim will speed up searching down below
x_mm = unique(avg_brain_mm_list(:,1));
y_mm = unique(avg_brain_mm_list(:,2));
z_mm = unique(avg_brain_mm_list(:,3));

%% put all the parcels into the space of the MNI brain
rois_tal = zeros(avg_brain_dims);
rois_anat = zeros(avg_brain_dims);

for ll = 1:length(labels)
   
    allcoords = find(parcels_nii.vol==labels(ll));
    [i,j,k]  = ind2sub(parcel_dims,allcoords);
    
    % list all the mm coordinates of voxels that have this label
    mm_coords = [xlist(i),ylist(j),zlist(k)];
    
    % convert these into anatomy space
%     mm_coords_anat = (mm_coords+ xfm_mat(:,4)')*xfm_mat(:,1:3) - xfm_mat(:,4)';
    mm_coords_anat = mm_coords*pinv(xfm_mat(:,1:3)) - xfm_mat(:,4)';
    for mm=1:length(mm_coords)
       
        % find all voxels that are close to this point (don't leave any
        % holes!)         
        xc = find(abs(x_mm-mm_coords(mm,1))<voxmm);
        yc = find(abs(y_mm-(-mm_coords(mm,3)))<voxmm);
        zc = find(abs(z_mm-mm_coords(mm,2))<voxmm);
      
        % label all of those voxels with the same label they had before
        [xc,yc,zc] = meshgrid(xc,yc,zc);
        xc = xc(:);yc = yc(:);zc = zc(:);
        for ii=1:numel(xc) 
            rois_tal(xc(ii),yc(ii),zc(ii)) = labels(ll);
        end
        
        % same for anatomical pts       
        xc = find(abs(x_mm-mm_coords_anat(mm,1))<voxmm);
        yc = find(abs(y_mm-(-mm_coords_anat(mm,3)))<voxmm);
        zc = find(abs(z_mm-mm_coords_anat(mm,2))<voxmm);
      
        % label all of those voxels with the same label they had before
        [xc,yc,zc] = meshgrid(xc,yc,zc);
        xc = xc(:);yc = yc(:);zc = zc(:);
        for ii=1:numel(xc) 
            rois_anat(xc(ii),yc(ii),zc(ii)) = labels(ll);
        end

    end
    
%     % convert this label into anat space
%     [pts_x,pts_y,pts_z] = ind2sub(size(rois_tal),find(rois_tal==labels(ll)));
%     
%     pts_anat = [pts_x, pts_y, pts_z]*xfm_mat;
end

%% save a new nifti with the MD definitions
new_nii = avg_brain;
new_nii.vol = rois_tal;
save_nifti(new_nii, parcels_xfm_nii);

new_nii = avg_brain;
new_nii.vol = rois_anat;
save_nifti(new_nii, parcels_anat_nii);

%% this part makes plots...a check to make sure the above worked. 
% mean_z = sum(sum(rois_tal==labels(9), 2),1);
% z2plot_mask = find(mean_z==max(mean_z));
% 
% figure;imagesc(rois_tal(:,:,z2plot_mask));colorbar();title('z-slice')
% figure;imagesc(avg_brain.vol(:,:,z2plot_mask));colorbar();title('z-slice')
% 
% mean_y = sum(sum(rois_tal==labels(9), 3),1);
% y2plot_mask = find(mean_y==max(mean_y));
% 
% figure;imagesc(squeeze(rois_tal(:,y2plot_mask,:)));colorbar();title('y-slice')
% figure;imagesc(squeeze(avg_brain.vol(:,y2plot_mask,:)));colorbar();title('y-slice')
% 
% mean_x = sum(sum(rois_tal==labels(9), 3),2);
% x2plot_mask = find(mean_x==max(mean_x),1);
% 
% figure;imagesc(squeeze(rois_tal(x2plot_mask,:,:)));colorbar();title('x-slice')
% figure;imagesc(squeeze(avg_brain.vol(x2plot_mask,:,:)));colorbar();title('x-slice')


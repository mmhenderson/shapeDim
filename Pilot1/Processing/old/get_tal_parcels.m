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

%% load the nifti with the MD ROIs masked out
% need to have a copy of this file in our current folder, this is the thing
% we downloaded from Fedorenko paper
parcels_nii = load_nifti(fullfile(exp_path, 'Processing', 'MDROI.nii.gz'));

% this is the name of the new file we will make here
parcels_xfm_nii = fullfile(exp_path, 'Processing', 'MDROI_XFM_TAL.nii.gz');

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
anat = load_nifti(tal_brain_orig_nii);
anat_dims = size(anat.vol);

% list the mm coordinates of every voxel in this array
% sform tells us a mapping from i,j,k, in the volume to x,y,z in mm
anat_sform = anat.sform;
[i,j,k] = meshgrid(1:anat_dims(1),1:anat_dims(2),1:anat_dims(3));
anat_mm_list = [i(:),j(:),k(:), ones(numel(i(:)),1)]*anat_sform';
anat_mm_list = anat_mm_list(:,1:3);

% taking unique in each dim will speed up searching down below
x_mm = unique(anat_mm_list(:,1));
y_mm = unique(anat_mm_list(:,2));
z_mm = unique(anat_mm_list(:,3));

%% put all the parcels into the space of the MNI brain
new_masked_vol = zeros(anat_dims);

for ll = 1:length(labels)
   
    allcoords = find(parcels_nii.vol==labels(ll));
    [i,j,k]  = ind2sub(parcel_dims,allcoords);
    
    % list all the mm coordinates of voxels that have this label
    mm_coords = [xlist(i),ylist(j),zlist(k)];

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
            new_masked_vol(xc(ii),yc(ii),zc(ii)) = labels(ll);
        end

    end
    
end

%% save a new nifti with the MD definitions
new_nii = anat;
new_nii.vol = new_masked_vol;
save_nifti(new_nii, parcels_xfm_nii);

%% this part makes plots...a check to make sure the above worked. 
mean_z = sum(sum(new_masked_vol==labels(9), 2),1);
z2plot_mask = find(mean_z==max(mean_z));

figure;imagesc(new_masked_vol(:,:,z2plot_mask));colorbar();title('z-slice')
figure;imagesc(anat.vol(:,:,z2plot_mask));colorbar();title('z-slice')

mean_y = sum(sum(new_masked_vol==labels(9), 3),1);
y2plot_mask = find(mean_y==max(mean_y));

figure;imagesc(squeeze(new_masked_vol(:,y2plot_mask,:)));colorbar();title('y-slice')
figure;imagesc(squeeze(anat.vol(:,y2plot_mask,:)));colorbar();title('y-slice')

mean_x = sum(sum(new_masked_vol==labels(9), 3),2);
x2plot_mask = find(mean_x==max(mean_x),1);

figure;imagesc(squeeze(new_masked_vol(x2plot_mask,:,:)));colorbar();title('x-slice')
figure;imagesc(squeeze(anat.vol(x2plot_mask,:,:)));colorbar();title('x-slice')


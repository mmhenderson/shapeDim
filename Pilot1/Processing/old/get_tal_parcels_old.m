%% load the anatomical definitions of a parcellation of the MD network 
% from Fedorenko, Duncan , & Kanwisher, 2013
% downloaded from imaging.mrc-cbu.cam.ac.uk/imaging/MDsystem

% raw data is a nifti in MNI space - convert to TAL space
% also save out MDROIs_tal.mat, that gets used by the next function

% MMH 6/30/17
%% 
clear
% mni2tal lives here
addpath('/usr/local/serenceslab/serenceslab_toolboxes/eeglab13_6_5b/plugins/dipfit2.3');


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

%% load the nifti 
parcels = load_nifti(fullfile(exp_path, 'Processing', 'MDROI.nii.gz'));


% range of the image in mm
% this is from the far edges of the voxels
% got these by loading in fslview
% x_range = [-78, 78];
x_range = [-79,79];
% y_range = [-112, 76];
y_range = [-113 , 77];
% z_range = [-50, 86];
z_range=  [-51, 87];

voxmm = 2;

% the coordinate of each pt is the center of this voxel in mm
xlist = (x_range(1)+1:voxmm:x_range(2));
ylist = (y_range(1)+1:voxmm:y_range(2));
zlist = (z_range(1)+1:voxmm:z_range(2));



labels = unique(parcels.vol(:));

% list the mm coordinates of every voxel in this array
parcels_sform = parcels.sform;
parcel_dims = size(parcels.vol);
[i,j,k] = meshgrid(1:parcel_dims(1),1:parcel_dims(2),1:parcel_dims(3));
parcel_mm_list = [i(:),j(:),k(:), ones(numel(i(:)),1)]*parcels_sform';
parcel_mm_list = round(parcel_mm_list(:,1:3));

if any(parcel_dims~=[length(xlist),length(ylist),length(zlist)])
    error('nifti image is wrong size')
end
%% convert the parcels into anatomical space, using my transformation matrix

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

%% load the anatomy (orig.nii.gz)

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

anat = load_nifti(reference_volume_nii);

% list the mm coordinates of every voxel in this array
anat_dims = size(anat.vol);
anat_sform = anat.sform;
[i,j,k] = meshgrid(1:anat_dims(1),1:anat_dims(2),1:anat_dims(3));
anat_mm_list = [i(:),j(:),k(:), ones(numel(i(:)),1)]*anat_sform';
anat_mm_list = round(anat_mm_list(:,1:3));

x_mm = unique(anat_mm_list(:,1));
y_mm = unique(anat_mm_list(:,2));
z_mm = unique(anat_mm_list(:,3));

%% generate a list of the coordinates for each parcel, in MNI coords (mm)
% then convert each mm coord to TAL
parcels_xfm_nii = fullfile(exp_path,'Processing','MDROI_XFM.nii.gz');
new_nii = anat;

% for ll=1:length(labels)
for ll=2
    
    allcoords = find(parcels.vol==labels(ll));
    [i,j,k]  = ind2sub(parcel_dims,allcoords);
    
    inds_mm = [xlist(i)',ylist(j)',zlist(k)'];
    
    tal_mm = mni2tal(inds_mm);
    
%     tal_mm = [i,j,k, ones(size(i))]*parcels_sform';
%     tal_mm = tal_mm(:,1:3);
    
    parcels(ll).inds_mm = inds_mm;
    parcels(ll).tal_mm = tal_mm;
    
    % now convert it to anatomical space, using tal transformation matrix
%     anat_mm = tal_mm*xfm_mat;
    anat_mm = tal_mm;
    anat_mm = anat_mm(:,1:3);
    
    parcels(ll).anat_mm = anat_mm;
    
    
    mask = zeros(anat_dims);
    
    for vv = 1:size(anat_mm,1)
       
%         vv
        xc = find(x_mm==round(anat_mm(vv,1)));
        yc = find(y_mm==round(anat_mm(vv,2)));
        zc = find(z_mm==round(anat_mm(vv,3)));
      
        mask(xc,yc,zc) = 1;
        
    end
    
    new_nii.vol = mask;
    
    
    parcels(ll).labstr = [num2str(ll) '-' num2str(labels(ll))];
    
end
%name of the file with all the parcels in TAL coords
newstruct = fullfile(exp_path,'Processing', 'MDROIs_tal.mat');

save(newstruct, 'parcels')

save_nifti(new_nii, parcels_xfm_nii);

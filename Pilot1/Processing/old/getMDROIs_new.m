% Transform all MD-network ROIs into subject functional space
% using FLIRT. 
% MD ROIs must be created in subject anatomical space before this step.

%%
clear
close all

sub = 'S01';
sub_doreti = 'CG';
exp_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';
doreti_path = '/usr/local/serenceslab/Doreti/';

unzip = 0;

%%

if ~exist([exp_path 'DataPreproc/' sub '/regheadercenterMOD.mat'],'file')
    error('need your registration matrix first!')
end

% Need a reference volume, this should be the first volume from your
% functional data (has to be session 1!!)
% reference_volume = [exp_path, 'DataPreproc/' sub '/MCTemplate01.nii.gz'];
reference_volume = [exp_path, 'DataPreproc/' sub '/MCTemplateXFM01.nii.gz'];

subpath = [exp_path, 'DataPreproc/', sub, '/']; % path to this subject's preprocessed data

labels = 1:28;
for ll = 1:length(labels)


    % get the brain mask file into nii format
    % file2align = [doreti_path 'ANAT/'  sub_doreti '/mri/brainmask.nii.gz'];
    file2align = fullfile(exp_path,'VOIs_MD',sub,sprintf('MDROI_label%d_REG2ANAT.nii.gz',ll));
    if ~exist(file2align,'file')
        error('need to align MNI brain to anat space first!')
    end

    file2save = fullfile(exp_path,'VOIs_MD',sub,sprintf('MDROI_label%d_REG2FUNC.nii.gz',ll));

    fprintf('loading your original ROI def from %s\n',file2align)

    fprintf('saving your transformed ROI def to %s\n',file2save)
    %%
    % make the transformation
    err =  unix(['flirt -in ', file2align ' -ref ', reference_volume,...
        ' -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', file2save]);
    if err
        error('your flirt command failed!')
    else
        disp('SUCCESS')
    end

%     old_nii = load_nifti(file2align);
    new_nii = load_nifti(file2save);
    new_nii.vol = double(new_nii.vol>0.5);
    save_nifti(new_nii,file2save);
end
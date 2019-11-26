% Make a brain mask file for your subject in the functional space
% corresponding to this experiment. This script finds a brain mask file
% that lives in Doreti folder in anatomical space, then uses your
% regheadercenterMOD.mat transformation matrix to map it into the
% functional space that your preprocessed niftis are in. Can then use it to
% mask out good voxels for further analyses!

% the resulting file is saved in DataPreproc/S##/BrainMask_REG2FUNC.nii

% always a good idea to check this mask in fsleyes just to make sure
% nothing is wrong!!

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

% get the brain mask file into nii format
% file2align = [doreti_path 'ANAT/'  sub_doreti '/mri/brainmask.nii.gz'];
file2align = fullfile(exp_path,'DataPreproc',sub,'MNI_BRAIN_MASK_XFMANAT.nii.gz');
if ~exist(file2align,'file')
    error('need to align MNI brain to anat space first!')
end

file2save = fullfile(exp_path,'DataPreproc',sub,'MNI_BRAIN_MASK_REG2FUNC.nii.gz');

fprintf('loading your original brain mask from %s\n',file2align)

fprintf('saving your transformed brain mask to %s\n',file2save)
%%
% make the transformation
err =  unix(['flirt -in ', file2align ' -ref ', reference_volume,...
    ' -applyxfm -init ', subpath, 'regheadercenterMOD.mat -out ', file2save]);
if err
    error('your flirt command failed!')
else
    disp('SUCCESS')
end


% and unzip it if you want to
if unzip
    fprintf('unzipping your transformed brain mask...\n')

    err = unix(['gzip -dk ' file2save]);
    if err
        error('your unzip command failed!')
    else
        disp('SUCCESS')
    end
end
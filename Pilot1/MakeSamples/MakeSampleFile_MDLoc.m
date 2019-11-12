% Make Sample File for shapeDim project - takes pre-processed niftis and
% extracts the voxels we care about.

clear
close all

% set inputs
FS_sbj = 'CG';
subnum = '01';
silent = 0;

% set paths (end all with filesep!)
experiment_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';

out_path = [experiment_path, 'Samples/'];
beh_path = [experiment_path, 'DataBehavior/'];
func_path = [experiment_path, 'DataPreproc/S', char(subnum), '/'];

if ~exist(out_path, 'dir'), mkdir(out_path); end
filename = fullfile(out_path, ['SampleFile_MDLoc_S', subnum]);

hemis = {'lh', 'rh'};

addpath('/usr/local/freesurfer6/matlab/')   % this folder has the load_nifti function we need

%% load the timing file (made in GetEventTiming.m)

fn = [out_path 'TimingFile_S' subnum '.mat'];
if ~exist(fn, 'file')
    error('need to make timing file first, run GetEventTiming.m')
end
fprintf('Loading event timing file\n')
load(fn)

%% Read VOIs
% Read the retino VOI files, which consist of a 3D volume with 1s marking
% the location of the ROI. We will store all the indices of each ROI, and
% then make sure none overlap. 

fprintf('Loading all ROIs\n')

% this will be a structure array [2 x nROIs] listing all my rois and their properties.
ROIs = struct('voxel_inds',[], 'name','');

for hh = 1:length(hemis)   
    
    % look at every VOI file in the folder.
    VOIfiles = dir([experiment_path, 'VOIs/S' char(subnum) '/', hemis{hh},'*.nii.gz']);

    % just retinotopic ones...these will have no extra extension at the
    % end.
    mdloc = contains({VOIfiles.name}, 'mdloc');
    
    VOIfiles(~mdloc) = [];
    
    for file_idx = 1:length(VOIfiles) 
        nifti = load_nifti([experiment_path, 'VOIs/S' char(subnum) '/', VOIfiles(file_idx).name]); 
        
        % extracting the indices from the binary mask with "find". now each
        % voxel gets a unique number that will stay with it throughout
        % preprocessing.
        voxel_inds_this_ROI = find(nifti.vol)';        
       
        %find where the filename starts with giving the file extension (eg
        %.nii.gz), and cut this filename extension off of it
        idx = min(strfind(VOIfiles(file_idx).name, '.nii')); 
        ROI_name = VOIfiles(file_idx).name(1:idx-1); 
        
        if hh==1
            % first hemisphere - every ROI gets a new column of struct array
            v_idx = file_idx;
        else
            % second hemisphere - check if there is a matching one in other
            % hemisphere, otherwise make a new column
            v_idx = find(contains({ROIs(1,:).name},ROI_name(3:end)));
            if isempty(v_idx)
               v_idx = size(ROIs,2)+1;
               % also initialize an empty name string for the other hemi
               ROIs(1,v_idx).name = '';
            end
        end
        
        % add information to the structure
        ROIs(hh,v_idx).voxel_inds = voxel_inds_this_ROI; 
        ROIs(hh,v_idx).name = ROI_name; 
%         ROIs(hh,v_idx).is_motor=contains(ROI_name,'DIGITLOC');

    end 
end 

nROIs = size(ROIs,2);

% make sure there's a string in every field, even if empty - this way we
% can use the contains function without error.
for hh=1:length(hemis)
    emptynames = find(cellfun(@isempty,{ROIs(hh,:).name}));
    for ii=1:length(emptynames)
        ROIs(hh,emptynames(ii)).name = '';
    end
end

%% Correct overlap
% Prune voxels which are shared between pairs of ROIs (interleaved):
for hh = 1:length(hemis)
    hemivox = [];
    for vv1 = 1:nROIs 
        %for the n visual areas in this hemisphere, compare that visual area against all other n-1 areas
        for vv2 = vv1+1:nROIs 
            
            % find the overlap
            overlapvox = intersect(ROIs(hh,vv1).voxel_inds, ROIs(hh,vv2).voxel_inds);
            if ~isempty(overlapvox)
                fprintf('detected %d voxels of overlap between %s and %s\n',numel(overlapvox),ROIs(hh,vv1).name,ROIs(hh,vv2).name);

                % split them up evenly.
                roi1_deletevox = overlapvox(1:2:end); %uneven voxels will deleted from roi1
                roi2_deletevox = overlapvox(2:2:end); %even voxels will be deleted from roi2                 
                ROIs(hh,vv1).voxel_inds(ismember(ROIs(hh,vv1).voxel_inds, roi1_deletevox)) = [];
                ROIs(hh,vv2).voxel_inds(ismember(ROIs(hh,vv2).voxel_inds, roi2_deletevox)) = [];
            end
        end       
    end
end

% Prune voxels which are shared between hemispheres (interleaved)
overlapvox = intersect([ROIs(1,:).voxel_inds], [ROIs(2,:).voxel_inds]);
hemi1_deletevox = overlapvox(1:2:end); 
hemi2_deletevox = overlapvox(2:2:end);
fprintf('correcting hemisphere overlap: found %d voxels\n',numel(overlapvox))
% delete these voxels from whichever ROI they belong to.
for vv = 1:nROIs

    todelete = intersect(ROIs(1,vv).voxel_inds, hemi1_deletevox);
    if ~isempty(todelete)
        ROIs(1,vv).voxel_inds(ismember(ROIs(1,vv).voxel_inds, todelete)) = [];
    end    
    
    todelete = intersect(ROIs(2,vv).voxel_inds, hemi2_deletevox);
    if ~isempty(todelete)
        ROIs(2,vv).voxel_inds(ismember(ROIs(2,vv).voxel_inds, todelete)) = [];
    end
end

% finally, unwrap all voxels that landed in any of these ROIs.
all_vox_concat = [ROIs.voxel_inds];
all_vox_concat = unique(all_vox_concat);

%% Get Sample Timecourse
% This will be a matrix of [nTimepoints x nVoxels]

RunsListMain = []; 

nTRs_main = 387 - 16;

fid = fopen([func_path, 'runs.list']);
line = fgetl(fid);
while ischar(line) && ~isempty(line)
    if strcmp(line(8:end),'stim')
        RunsListMain = [RunsListMain; line(1:6)];
    end
   
    line = fgetl(fid);    
end
fclose(fid);

fprintf('loading main task niftis:\n')
% Make samplesMain
tmp=NaN(nTRs_main,length(all_vox_concat),length(RunsListMain));
for run = 1:length(RunsListMain) %for each functional run load the corresponding nifti
    fprintf('    %d of %d\n',run,length(RunsListMain));
    nifti_run = load_nifti([func_path, RunsListMain(run,1:2), '_REG_MC_DET_', RunsListMain(run,4:6), '.nii.gz']);
    reshaped_volume = reshape(nifti_run.vol, [prod(nifti_run.dim(2:4)) nifti_run.dim(5)])'; %Reshape so that we have one index for each voxel (corresponding to the indices in our VOI file), with voxels as columns
    tmp(:,:,run)=reshaped_volume(:,all_vox_concat);
end

samplesMain = [];
for run = 1:length(RunsListMain)
    samplesMain=[samplesMain;tmp(:,:,run)];
end

clear tmp



%% Save Sample File


fprintf('saving all data to %s\n',filename);
save(filename, 'ROIs','all_vox_concat','main','samplesMain','-v7.3');


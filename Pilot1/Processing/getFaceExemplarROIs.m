% Define ROIs that show selectivity for particular exemplars using the
% neutral task (result of a one-way anova)

% we'll use the freesurfer automatic cortical parcellation to define rough anatomical parcels, then these
% get thresholded with the face localizer and saved as final ROIs.

% need to run the scripts in AnalyzeNeutralTask first.

%%
clear 
close all

subnum = '01';
subinit = 'AV';

pthresh = 0.05;

% these two contrasts have to correspond to copes 1 and 2!
% contrasts = {'left>right','right>left'};
hemis = {'lh','rh'};
nHemis = length(hemis);

retfolder = '/usr/local/serenceslab/Doreti/';
subretfolder = [retfolder 'ANAT/' subinit '/'];
exp_path = '/mnt/neurocube/local/serenceslab/maggie/faceDim/pilot4/';
func_path = [exp_path, 'DataPreproc/S', char(subnum), '/'];

% load the annotation file so that we can get a list of all the names (this
% can really be any subject or hemisphere because they're all the same
% names)
% next we will make a volume for every single
% cortical area, then see how many intersect with the face>place localizer. 
[inds,labs,info]= read_annotation('/usr/local/serenceslab/Doreti/ANAT/BX/label/lh.aparc.annot');
allnames = info.struct_names(2:end);
ROIs = allnames;
ROIs_4unix = ROIs;
ROIs_4saving = ROIs;
for rr=1:length(ROIs)
    if contains(ROIs{rr},'&')
        % put an extra slash in here so that unix won't be confused by the & signs
        newname = ROIs{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '\' newname(andind:end)];
        ROIs_4unix{rr} = newname;
        newname = ROIs{rr};
        andind = find(newname=='&');
        newname = [newname(1:andind-1) '_and_' newname(andind+1:end)];
        ROIs_4saving{rr} = newname;
    end
end
clear info

% this extension gets appended to the final nifti files for each volume,
% indicating they were thresholded with this localizer.
ROI_ext  = '_EXLOC';

% this first dir is a place to store the un-thresholded VOIs after they are
% made, the second dir is our usual VOI folder
VOIdir1 = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
VOIdir2 = [exp_path, 'VOIs/S' char(subnum) '/'];

minVox = 50;

currentdir = pwd;

% set this to 1 if you want to make all the volumes in functional space.
% takes a while!
doVolumes = 1;

%% first go into the Doreti folder, and make labels based on the anatomy

if doVolumes

    % check to see if labels are already made or not
    makelabs = true;
    if exist([subretfolder, 'label/parcels_all'], 'dir')
        if ~isempty(dir([subretfolder, 'label/parcels_all/*.label'])) % if there are labels in here already
            rsp = input('Labels already exist in this subject''s parcels_all folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makelabs = false;
            end
        end
    else
        mkdir([subretfolder, 'label/parcels_all']);
    end

    if makelabs
        cd(retfolder)
        for hh=1:nHemis
 %             unix(['mri_annotation2label --annotation aparc.a2009s --subject ' subinit ' --hemi ' hemis{hh} ' --outdir ' subretfolder 'label/parcels_all/']); 
            unix(['mri_annotation2label --annotation aparc --subject ' subinit ' --hemi ' hemis{hh} ' --outdir ' subretfolder 'label/parcels_all/']); 

        end
        cd(currentdir)
    end
    %% now turn them into VOI masks in functional space

    % here we are making a new folder to hold ROIs that were generated using
    % fsl parcellation scheme. These have to be processed a bit more before
    % they can be put into our normal pipeline. 

    % check to see if masks are already made or not
    makemasks = true;
    if exist([exp_path, 'VOIs_parcel/S', char(subnum)], 'dir')
        if ~isempty(dir([exp_path, 'VOIs_parcel/S', char(subnum), '/*h_*nii.gz'])) % if there are visual masks found
            rsp = input('Visual area masks already exist in this subject''s VOIs_parcel folder. Overwrite? (y/n) ', 's');
            if strcmp(rsp, 'n')
                makemasks = false;
            end
        end
    else
        mkdir([exp_path, 'VOIs_parcel/S', char(subnum)]);
    end

    % make Visual Area masks
    if makemasks
        delete([exp_path, 'VOIs_parcel/S' char(subnum) '/*h_*nii.gz']); %delete retino content from VOI folder 
    %     ROIs_orig = {'S_precentral-sup-part','S_precentral-inf-part','S_central','S_postcentral',...
    %         'G_precentral','G_postcentral','G_parietal_sup',...
    %         'G&S_subcentral','G_pariet_inf-Supramar','G&S_paracentral','G_front_sup'};

        func_template = [func_path, 'MCTemplateXFM01.nii.gz'];
        regfile = [func_path, 'Func2Anat_auto.dat'];
        outdir = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
        if ~isdir(outdir)
            mkdir(outdir)
        end
        if ~exist(regfile,'file')
            error('registration to anat not made yet! you need to run getVisualROIs.m first')
        end
        % below is in replacement of label2VOI function
        for hemi_idx = 1:length(hemis)
            for roi_idx = 1:length(ROIs_4unix)
                labelfile = [subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', ROIs_4unix{roi_idx}, '.label'];


                if exist([subretfolder, '/label/parcels_all/', hemis{hemi_idx}, '.', ROIs{roi_idx}, '.label'], 'file')

                    unixcmd = ['mri_label2vol --label ', labelfile, ' --temp ', func_template, ' --subject ', subinit, ' --hemi ', hemis{hemi_idx}, ...
                        ' --o ', outdir, hemis{hemi_idx}, '_', ROIs_4unix{roi_idx}, '.nii.gz --proj frac 0 1 0.1 --fillthresh 0.3 --reg ', regfile];
                    [s,~] = unix(unixcmd, '-echo');
                    if s==1
                        out(['\nERROR: An error occured while processing label ', labelfile, '. Exiting.\n']);
                        status = 1;
                        return
                    end                
                end            
            end %roi_idx   
        end
    end    
end

%% now load the localizer, and mask out the parts we need

% re-organize them into a table [nRuns x nVox] including only voxels in
% vInd list
pstats_nii = load_nifti([exp_path 'AnalyzeNeutralTask/S' subnum '/ExemplarSelectivity_p.nii']);
pstats = pstats_nii.vol;
% dims = size(pstats);
% pstats = reshape(pstats,prod(dims(1:3)),dims(4))';

% this is a list of all the voxels that passed threshold 
% voxels_thresholded = find(pstats<pthresh);
vals2use = pstats<pthresh;

fprintf('over whole brain, %d voxels passed the F-test threshold on full training data set\n',sum(vals2use(:))), 
% 
% fstats_nii = load_nifti([experiment_path 'AnalyzeNeutralTask/S' subnum '/ExemplarSelectivity_F.nii']);
% fstats = fstats_nii.vol(:);
% dims = size(fstats);
% fstats = reshape(fstats,prod(dims(1:3)),dims(4))';

%% define my final mask

% vals2use = voxels_thresholded;

%% create each VOI file

for hh=1:length(hemis)
     
    for vv=1:length(ROIs)

        if exist([subretfolder, '/label/parcels_all/', hemis{hh}, '.', ROIs{vv}, '.label'], 'file')

            fullVOI = [VOIdir1 hemis{hh} '_' ROIs_4unix{vv} '.nii.gz'];
            thisVOI = load_nifti(fullVOI);

            fullVol = thisVOI.vol;

            maskedVol = fullVol & vals2use;
    %         fprintf('%s %s %s: found %d voxels for %s\n',subinit,hemis{hh}, ROIs{vv},sum(maskedVol(:)),contrasts{cc});

            if sum(maskedVol(:))>minVox
                fprintf('%s %s %s: found %d voxels with exemplar selecivity in neutral task\n',subinit,hemis{hh}, ROIs{vv},sum(maskedVol(:)));
                newVOI = thisVOI;
                newVOI.vol = maskedVol;
                savename = [VOIdir2 hemis{hh} '_' ROIs_4saving{vv} ROI_ext '.nii.gz'];
                err = save_nifti(newVOI,savename);
                if err~=0
                    error('error saving new VOI file')
                end
            end
        end

    end
        
end
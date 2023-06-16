function [] =  make_anat_VOIs(sub2do)

% Define ROIs that show selectivity for object>scrambled with the category localizer.

% we'll use the freesurfer automatic cortical parcellation to define rough anatomical parcels, then these
% get thresholded with the face localizer and saved as final ROIs.

% Need to run GLM using feat to create the t-map first.

if nargin<1
    sub2do=[1,2,3,4,5,6,7];
end

subinit_big = {'CP','BX','BR','CA','CG','CT','CU'};
subnum_big = [1,2,3,4,5,6,7];

for ss = sub2do

    % set inputs
    subinit = subinit_big{ss};
    subnum = sprintf('%02d',subnum_big(ss));

    % these two contrasts have to correspond to copes 1 and 2!
    % contrasts = {'left>right','right>left'};
    hemis = {'lh','rh'};
    nHemis = length(hemis);

    retfolder = '/usr/local/serenceslab/Doreti/';
    subretfolder = [retfolder 'ANAT/' subinit '/'];
    exp_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/';
    func_path = [exp_path, 'DataPreproc/S', char(subnum), '/'];

    % load the annotation file so that we can get a list of all the names (this
    % can really be any subject or hemisphere because they're all the same
    % names)
    % next we will make a volume for every single
    % cortical area, then see how many intersect with the face>place localizer. 
    [inds,labs,info]= read_annotation(['/usr/local/serenceslab/Doreti/ANAT/', ...
        subinit '/label/lh.aparc.annot']);
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


    % cope 1 is face>place, cope 7 is object>scrambled
    % cope2use = [1, 5];
    % copenames = {'facecatloc','objcatloc'};
    cope2use = [5];
    copenames = {'objcatloc'};

    % this first dir is a place to store the un-thresholded VOIs after they are
    % made, the second dir is our usual VOI folder
    VOIdir1 = [exp_path, 'VOIs_parcel/S' char(subnum) '/'];
    VOIdir2 = [exp_path, 'VOIs/S' char(subnum) '/'];

    minVox = 50;

    currentdir = pwd;


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
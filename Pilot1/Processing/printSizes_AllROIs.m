% how many voxels are in each of my ROIs?

%%
clear 
close all

root = '/usr/local/serenceslab/maggie/faceDim/pilot4/';

sublist = {'01'};

addpath('/usr/local/freesurfer6/matlab/') % add my load_nifti function

%% retinotopic ROIs first

my_areas_orig = {'V1v','V2v','V3v','V1d','V2d','V3d','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};

my_areas = my_areas_orig;

hemis = {'lh','rh'};

all_sizes = zeros(length(my_areas), length(hemis),length(sublist));

col_names = [];

for ss=1:length(sublist)
    
    subnum = sublist{ss};
    
    col_names = [col_names, {[hemis{1} '_' subnum], [hemis{2} '_' subnum]}];    

    for vv=1:length(my_areas)
        
        for hh=1:length(hemis)
            
            fullname = [root, 'VOIs/S' char(subnum) '/', hemis{hh} '_' my_areas{vv} '.nii.gz'];
    
            if ~exist(fullname,'file')
                fprintf('%s does not exist\n',fullname)
                all_sizes(vv,hh,ss) = 0;
                continue
            end
%             fprintf('loading %s\n',fullname);
            
            nifti = load_nifti(fullname); %load the nifti with the binary mask
            
            ROI_voxels = find(nifti.vol); %and get the indices of voxels that are in this visual area

            all_sizes(vv,hh,ss) = numel(ROI_voxels);

        end
    end
    
    
end

fprintf('Retinotopic ROI sizes:\n\n');
    
tab = array2table(reshape(all_sizes, length(my_areas), length(hemis)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_orig);

disp(tab)

%% Face>scrambled selective areas:

hemis = {'lh','rh'};

all_sizes = [];

col_names = [];

my_areas = [];

for ss=1:length(sublist)
    
    subnum = sublist{ss};
    
    col_names = [col_names, {[hemis{1} '_' subnum], [hemis{2} '_' subnum]}];    
    
    for hh = 1:length(hemis)   

        % look at every VOI file in the folder.
        VOIfiles = dir([root, 'VOIs/S' char(subnum) '/', hemis{hh},'*CATLOC.nii.gz']);

        if isempty(VOIfiles)
            all_sizes(1,hh,ss) =0;
        end
        for file_idx = 1:length(VOIfiles) 
            nifti = load_nifti([root, 'VOIs/S' char(subnum) '/', VOIfiles(file_idx).name]); 

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
                v_idx = find(contains(my_areas,ROI_name(4:end)));
                if isempty(v_idx)
                   v_idx = size(my_areas,2)+1;                  
                end
                
            end
            
            my_areas{v_idx} = ROI_name(4:end);
            all_sizes(v_idx,hh,ss) = numel(voxel_inds_this_ROI);
             
            
        end
    end 
end

fprintf('Face>Place ROI sizes:\n\n');
    
tab = array2table(reshape(all_sizes, length(my_areas), length(hemis)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas);

disp(tab)

%% Face exemplar selective areas:

hemis = {'lh','rh'};

all_sizes = [];

col_names = [];

my_areas = [];

for ss=1:length(sublist)
    
    subnum = sublist{ss};
    
    col_names = [col_names, {[hemis{1} '_' subnum], [hemis{2} '_' subnum]}];    
    
    for hh = 1:length(hemis)   

        % look at every VOI file in the folder.
        VOIfiles = dir([root, 'VOIs/S' char(subnum) '/', hemis{hh},'*EXLOC.nii.gz']);

        for file_idx = 1:length(VOIfiles) 
            nifti = load_nifti([root, 'VOIs/S' char(subnum) '/', VOIfiles(file_idx).name]); 

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
                v_idx = find(contains(my_areas,ROI_name(4:end)));
                if isempty(v_idx)
                   v_idx = size(my_areas,2)+1;                  
                end
                
            end
            
            my_areas{v_idx} = ROI_name(4:end);
            all_sizes(v_idx,hh,ss) = numel(voxel_inds_this_ROI);
             
            
        end
    end 
end

fprintf('Face Exemplar Selective (neutral task) ROI sizes:\n\n');
    
tab = array2table(reshape(all_sizes, length(my_areas), length(hemis)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas);

disp(tab)
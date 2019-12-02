% print out a table describing the sizes of all ROIs in the multiple-demand
% network (made using getMDROIs.m)

%%
clear 
close all

root = '/usr/local/serenceslab/maggie/shapeDim/Pilot1/';

% these are the subs that i have MD definitions for
sublist = {'01'};
my_areas_orig = {'IFS', 'AI-FO', 'iPCS', 'sPCS','poCS'};
my_areas = my_areas_orig;


for aa=1:length(my_areas)
   my_areas{aa}(my_areas{aa}==' ') = '_'; 
   my_areas{aa} = [my_areas{aa} '_mdloc']; 
end

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
            fprintf('loading %s\n',fullname);
            
            nifti = load_nifti(fullname); %load the nifti with the binary mask
            
            ROI_voxels = find(nifti.vol); %and get the indices of voxels that are in this visual area
            all_sizes(vv,hh,ss) = numel(ROI_voxels);

        end
    end
    
    
end

fprintf('Final ROI sizes:\n\n');
    
tab = array2table(reshape(all_sizes, length(my_areas), length(hemis)*length(sublist)), 'VariableNames',col_names, 'RowNames',my_areas_orig);

disp(tab)
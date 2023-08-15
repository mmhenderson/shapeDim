%% Generate shape space for ShapeDim fMRI pilot 1
% 4x4 grid only

% This version allows for background color to be chosen - good for the
% scanner because black shows dust from projector
%%
clear

close all;

addpath(genpath(fullfile(pwd,'ODBstims')));

kinds = {'balls','snm','fetus','cater'};

my_back_color = 0.3;
my_shape_color = 0.9;

%% loop over stim sets
for kk=1:length(kinds)
   
    savedir = fullfile(pwd,sprintf('AmpGrid%d_grey',kk));
    if ~isdir(savedir)
        mkdir(savedir);
    end
    kind = kinds{kk};
    mysize = 3000;
    stimtype = '';
    start = 0.2;    % min value along each axis
    stop = 4.8;
    
    % how many steps in the main grid (without the variable points?)
    nsteps_main = 4;
    main_pts = round(linspace(start,stop, nsteps_main),1);
    all_pts = main_pts; 
    [gridx,gridy] = meshgrid(all_pts,all_pts);
    all_grid_points = [gridx(:),gridy(:)];

    %% now actually craete and save the shapes
    for pp=1:size(all_grid_points,1)

        X=RFCcalc(ODBparams(kind),mysize,stimtype,all_grid_points(pp,:));

        image = X{1};
        image(image==0) = my_back_color;
        image(image==1) = my_shape_color;

        savepath = fullfile(savedir, sprintf('Shape_%.2f_%.2f.png',all_grid_points(pp,1),all_grid_points(pp,2)));
        imwrite(image,savepath);

        fprintf('saved image to %s\n',savepath);

    end

end
%% Generate full shape spaces for ShapeDim experiment

% this script will save out a larger grid of shapes, to allow for
% parametric manipulation of the difficulty. 
% saves about 2000 images (47^2) per group

% This version allows for background color to be chosen - good for the
% scanner because black shows dust from projector (trying to match these
% parameters closely for behav experiment)

% This set is to be used for behavioral pilot V2, stim set 3
% (starting data collection week of 10/07/2019)

%%
clear

close all;

addpath(genpath(fullfile(pwd,'ODBstims')));

kinds = {'balls','snm','fetus','cater'};

my_back_color = 0.3;
my_shape_color = 0.9;

%% loop over stim sets
for kk=[1:4]

    savedir = fullfile(pwd,sprintf('AmpGrid%d_full_grey_small',kk));
    if ~isdir(savedir)
        mkdir(savedir);
    end
    kind = kinds{kk};
    mysize = 224;
    stimtype = '';
    start = 0;    % min value along each axis
    stop = 5;
    step = 0.1;
    
    all_pts = start:step:stop;  
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
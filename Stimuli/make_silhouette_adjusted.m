%% overlay images and make a silhouette image
% save the image into the correct folder
%%
close all;
root = pwd;
kinds = {'balls','snm','fetus','cater'};

nsteps = 4;
start = 0.1;
stop = 4.9;

% here the "coordinates" for axes 1 and 2 are same spacing - the scaling of
% axis 1 has already been accounted for when making these images in shape
% space (make_grid_full_grey_adjusted.m), so all we want to do here is load
% and make the silhouette from these images.
all_pts_dim2 = linspace(start,stop,nsteps);  
all_pts_dim1 = all_pts_dim2;
[gridx,gridy] = meshgrid(all_pts_dim1,all_pts_dim2);
points = round([gridx(:),gridy(:)],1);

my_back_color = 0.3;
my_shape_color = 0.9;

% loop over stim sets
for kk=[3]

    folder = fullfile(root, sprintf('AmpGrid%d_adj_full_grey',kk));
   
    %% loop over grid points
    for xx=1:size(points,1)

        im_fn = fullfile(folder, sprintf('Shape_%.2f_%.2f.png',points(xx,1),points(xx,2)));
        im = imread(im_fn);
        if xx==1
             im_stack = zeros(size(im,1),size(im,2),size(points,1));
        end
        im_stack(:,:,xx) = im;
       
    end
    
    image = double(any(im_stack~=77, 3));
    image(image==1) = my_shape_color;
    image(image==0) = my_back_color;
    
    savepath = fullfile(folder, 'Silhouette_any.png');
    imwrite(image,savepath);

    fprintf('saved image to %s\n',savepath);
    
%     figure;hold all;
%     imshow(image);
% %     set(gca,'YDir','rev')
%     title(sprintf('set %d',kk))
end

%% overlay images and make a silhouette image
% save the image into the correct folder
%%
close all;
root = pwd;
kinds = {'balls','snm','fetus','cater'};

nsteps = 4;
start = 0.1;
stop = 4.9;
[gridx,gridy] = meshgrid(linspace(start,stop,nsteps),linspace(start,stop,nsteps));
points = round([gridx(:),gridy(:)],1);

my_back_color = 0.3;
my_shape_color = 0.9;

% loop over stim sets
for kk=[1:4]

    folder = fullfile(root, sprintf('AmpGrid%d_full_grey',kk));
   
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

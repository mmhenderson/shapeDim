%% Draw a grid of thumbnail images at various points in the shape space
% save out an image Grid.png, inside the relevant AmpGrid folder

%%
close all;

kinds = {'balls','snm','fetus','cater'};

% loop over stim sets
for kk=3

    savedir = fullfile(pwd,sprintf('AmpGrid%d',kk));
    kind = kinds{kk};
    mysize = 500;
    stimtype = '';
    nsteps = 4;
    start = 0.2;
    stop = 4.8;
    [gridx,gridy] = meshgrid(linspace(start,stop,nsteps),linspace(start,stop,nsteps));
    points = [gridx(:),gridy(:)];
    X=RFCcalc(ODBparams(kind),mysize,stimtype,points);

    figure('Position',[0,0,1000,1000]);
    hold all;axis off;axis equal
    % ax = [];
    axspace_ratio = 4;    % ratio of space to axes
    axsize = axspace_ratio/(axspace_ratio*nsteps+nsteps+1);
    assert(axsize*nsteps+axsize/axspace_ratio*(nsteps+1) ==1)
    axspace= axsize/axspace_ratio;

    axpos = axspace:axsize+axspace:1-axspace;
    [axpos_x,axpos_y] = meshgrid(axpos, axpos);

    axpos_full = [axpos_x(:), axpos_y(:), repmat(axsize,numel(axpos_x), 2)];

    %% loop over grid points
    for pp=1:size(points,1)

        ax = axes('Position',axpos_full(pp,:));hold all;

        thisim = X{pp};
        thisim_tmp = thisim;

        imshow(thisim);
        title(sprintf('[%.1f, %.1f]',points(pp,1),points(pp,2)));

    end
    %% end loop and save the image itself
    mi = 0.1;ma=1;cent = 0.51;
   
    annotation('textarrow',[cent,cent],[mi,ma],'String',sprintf('boundary %d',1),'Color','k','LineWidth',3)
    
    annotation('textarrow',[mi,ma],[cent,cent],'String',sprintf('boundary %d',2),'Color','k','LineWidth',3)
    
    set(gcf,'Color','w');
    cdata = print('-RGBImage','-r300');
    savepath = fullfile(savedir, sprintf('Grid_boundary%d.png',bb));
    imwrite(cdata, savepath);
    fprintf('saved image to %s\n',savepath);
    % sgtitle(sprintf('shape space %d',kk))

end
%% Measure the pixelwise similarity between images in each grid

%%
close all;

kinds = {'balls','snm','fetus','cater'};

euc_dist = zeros(length(kinds),16,16);
correl = ones(length(kinds),16,16);

all_ims = [];
kind_list = [];
coord_list = [];
% loop over stim sets
for kk=1:4

%     savedir = fullfile(pwd,sprintf('AmpGrid%d',kk));
    kind = kinds{kk};
    mysize = 100;
    stimtype = '';
    nsteps = 4;
    start = 0.2;
    stop = 4.8;
    [gridx,gridy] = meshgrid(linspace(start,stop,nsteps),linspace(start,stop,nsteps));
    points = [gridx(:),gridy(:)];
    X=RFCcalc(ODBparams(kind),mysize,stimtype,points);

    %% loop over grid points
    for pp1=1:size(points,1)
        
        all_ims = [all_ims; X{pp1}(:)'];
        kind_list = [kind_list; kk];
        coord_list = [coord_list; points(pp1,:)];
        
        for pp2=pp1+1:size(points,1)
            
            im1 = X{pp1}(:);
            im2 = X{pp2}(:);
            
            dist = sqrt(sum((im1-im2).^2));
            
            euc_dist(kk,pp1,pp2) = dist;
            euc_dist(kk,pp2,pp1) = dist;
            
            c = corrcoef(im1,im2);
            correl(kk,pp1,pp2) = c(2,1);
            correl(kk,pp2,pp1) = c(2,1);
        end

    end
end

%% plot dissimilarity matrix
cax = [];
for kk=1:4
    figure;hold all;
    cax = [cax, gca];
    
    imagesc(squeeze(euc_dist(kk,:,:)));
    colorbar();
    title(sprintf('Euc. distance, image set %d',kk));
end
match_clim(cax);
%% plot correlation matrix 
cax = [];
for kk=1:4
    figure;hold all;
    cax = [cax, gca];
    
    imagesc(squeeze(correl(kk,:,:)));colorbar();
    title(sprintf('Correlation, image set %d',kk));
end
match_clim(cax);


%% MDS

for kk=1:4
   
    D = squeeze(euc_dist(kk,:,:));
    
    reduced_rep = mdscale(D,2);
    
    figure;hold all;
    legendlabs = {};
    cols = parula(size(reduced_rep,1))
    for pp=1:size(reduced_rep,1)
        
        plot(reduced_rep(pp,1),reduced_rep(pp,2),'o','MarkerFaceColor',cols(pp,:),'MarkerEdgeColor',cols(pp,:))
        legendlabs{pp} = sprintf('%.1f, %.1f',points(pp,1),points(pp,2));
    end
    legend(legendlabs,'Location','EastOutside')
    title(sprintf('MDS on images from set %d\n(%d pix per image)',kk,numel(im1)));
    xlabel('MDS axis 1');
    ylabel('MDS axis 2');
    
end


%% PCA

[a,b] = pca(all_ims);

%%
figure;hold all;
for kk=1:4
    inds = kind_list==kk;
    plot(b(inds,1),b(inds,2),'o')
end
legend(kinds)
title('PCA on all images (16 from each set)');
xlabel('PC 1');
ylabel('PC 2');

%% 
cols = hsv(16);
for kk=1:4
    figure;hold all;
    legendlabs = {};
    for pp=1:size(points,1)
        inds = kind_list==kk & ismember(coord_list,points(pp,:),'rows');
        plot(b(inds,1),b(inds,2),'o','MarkerFaceColor',cols(pp,:))
        legendlabs{pp} = sprintf('%.1f, %.1f',points(pp,1),points(pp,2));
    end
    legend(legendlabs,'Location','EastOutside')
    title(sprintf('PCA on images from set %d\n(%d pix per image)',kk,numel(im1)));
    xlabel('PC 1');
    ylabel('PC 2');
end

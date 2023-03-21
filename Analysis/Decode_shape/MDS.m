% measure overall decodability of shapes

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1)-1);
% savepath = fullfile(root, 'Analysis','Decoding_results','PairwiseAcc_all.mat');

nStimsMain = 16;
nTasks = 3;
task_names = {'Linear (1+4 vs. 2+3)','Linear (1+2 vs. 3+4)','Checker (1+3 vs. 2+4)','Repeat detection'};

pairs = combnk(1:nStimsMain, 2);
nPairs = size(pairs,1);
     
pairwise_dist = zeros(size(pairs,1),1);

center = 2.5;
[xpoints,ypoints] = meshgrid(round(linspace(0.1, 4.9, 4),1), round(linspace(0.1, 4.9, 4),1));
xpoints = xpoints(:);
ypoints = ypoints(:);
points = [xpoints, ypoints];

quadrant = zeros(size(points,1),1);
quadrant(points(:,1)>center & points(:,2)>center) = 1;
quadrant(points(:,1)<center & points(:,2)>center) = 2;
quadrant(points(:,1)<center & points(:,2)<center) = 3;
quadrant(points(:,1)>center & points(:,2)<center) = 4;

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        euc_dist = zeros(nSubj, nVOIs, nTasks+1, nStimsMain, nStimsMain);
    end
    

    for vv=1:nVOIs

        for tt=1:nTasks
            
            fprintf('processing %s, task %d\n',ROI_names{vv}, tt);
            
            inds2use = mainSig(vv).BoundLabels==tt & mainSig(vv).IsMainLabels==1;
            
            pointlabs = mainSig(1).PointLabels(inds2use,:);
            [pts, ~, stim_inds] = unique(pointlabs,'rows');
            assert(all(all(pts==points)));

            dat = mainSig(vv).dat_avg(inds2use,:);
            dat = dat - repmat(mean(dat,2),1,size(dat,2));
            
            for ss1 = 1:nStimsMain

                for ss2 = ss1+1:nStimsMain

                    dat1 = dat(stim_inds==ss1,:);
                    dat2 = dat(stim_inds==ss2,:);

                    dist = get_normEucDist(dat1,dat2);

                    euc_dist(ss,vv,tt,ss1,ss2) = dist;
                    euc_dist(ss,vv,tt,ss2,ss1) = dist;
                    
                end
            end 
            
        end
    end
    
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    tt=nTasks+1;
    
    for vv=1:nVOIs

        fprintf('processing %s, repeat task\n',ROI_names{vv});

        inds2use = repSig(vv).IsMainLabels==1;

        pointlabs = repSig(1).PointLabels(inds2use,:);
        [pts, ~, stim_inds] = unique(pointlabs,'rows');
        assert(all(all(pts==points)));

        dat = repSig(vv).dat_avg(inds2use,:);
        dat = dat - repmat(mean(dat,2),1,size(dat,2));

        for ss1 = 1:nStimsMain

            for ss2 = ss1+1:nStimsMain

                dat1 = dat(stim_inds==ss1,:);
                dat2 = dat(stim_inds==ss2,:);

                dist = get_normEucDist(dat1,dat2);

                euc_dist(ss,vv,tt,ss1,ss2) = dist;
                euc_dist(ss,vv,tt,ss2,ss1) = dist;

            end
        end 

    
    end
end

%%
close all
ss=1;
% cols = plasma(nStimsMain);
cols = plasma(5);
cols = cols(1:4,:);

v2plot = 12:17;
for vv=1:numel(v2plot)
    
    figure;hold all;
    for tt=1:nTasks+1
        
        subplot(2,2,tt);hold all;
    

        D = squeeze(euc_dist(ss,v2plot(vv),tt,:,:));

        reduced_rep = mdscale(D,2);


        legendlabs = {};
        for qq=1:4

            plot(reduced_rep(quadrant==qq,1),reduced_rep(quadrant==qq,2),'o','MarkerFaceColor',cols(qq,:),'MarkerEdgeColor','None')
            legendlabs{qq} = sprintf('Q%d',qq);
        end
        legend(legendlabs,'Location','EastOutside')
        title(task_names{tt});
        xlabel('MDS axis 1');
        ylabel('MDS axis 2');
        axis('square')
        
    end
    suptitle(sprintf('MDS: %s\n',ROI_names{v2plot(vv)}));
    set(gcf,'Position',[200,200,1000,1000])
    set(gcf,'Color','w');
end
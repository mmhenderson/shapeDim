% measure overall decodability of shapes

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1));
savepath = fullfile(root, 'Analysis','Decoding_results','PairwiseAcc_all.mat');

%% setup the grid
nTasks = 3;
task_names = {'Linear (1)','Linear (2)','Checker','Repeat Detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);

nStimsMain = 16;
nRunsTotal = nTasks*2*2*3;

pairs = combnk(1:nStimsMain, 2);
nPairs = size(pairs,1);
     
pairwise_dist = zeros(size(pairs,1),1);

center_pos = [2.5, 2.5];
[xpoints,ypoints] = meshgrid(round(linspace(0.1, 4.9, 4),1), round(linspace(0.1, 4.9, 4),1));
xpoints = xpoints(:);
ypoints = ypoints(:);
points = [xpoints, ypoints];

center=center_pos(1);
quadrant = zeros(size(points,1),1);
quadrant(points(:,1)>center & points(:,2)>center) = 1;
quadrant(points(:,1)<center & points(:,2)>center) = 2;
quadrant(points(:,1)<center & points(:,2)<center) = 3;
quadrant(points(:,1)>center & points(:,2)<center) = 4;


which_dim = zeros(size(points,1),1);
across_bound = zeros(size(points,1),2);
for pp=1:size(pairs,1)
    pairwise_dist(pp) = sqrt(sum((points(pairs(pp,1),:) - points(pairs(pp,2),:)).^2));
    
    pt1 = points(pairs(pp,1),:);
    pt2 = points(pairs(pp,2),:);
    % do the points share the same y coordinate? if so, the discrimination
    % between them is along the x-coordinate.
    if pt1(2)==pt2(2)
        which_dim(pp) = 1;
    elseif pt1(1)==pt2(1)
        which_dim(pp) = 2;
    end
    % note that they can't ever have both coords shared...because they're all
    % different points. but they could have neither shared, which would
    % leave which_dim=0
    
    % do the points cross bounds 1 and 2? mark each column separately
    for dd =1:2
        across_bound(pp,dd) = (pt1(dd)>center_pos(dd)) + (pt2(dd)>center_pos(dd))==1;
    end
end

pairwise_dist_round = round(pairwise_dist*2)./2;

[unvals, ia,ib] = unique(pairwise_dist_round);

numeach = sum(repmat(ib, 1, numel(unvals))==repmat((1:numel(unvals)), numel(ib,1), 1), 1);


%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks+1, nPairs);
        acc_matrix = zeros(nSubj,nVOIs,nTasks+1,nStimsMain,nStimsMain);

    end

    for vv=1:nVOIs

        inds2use = mainSig(vv).IsMainLabels==1;
        assert(sum(inds2use)==32*nRunsTotal);
        
        dat = mainSig(vv).dat_avg(inds2use,:);

        pointlabs = mainSig(vv).PointLabels(inds2use,:);
        [pts, ~, stim_inds] = unique(pointlabs,'rows');
        assert(all(all(pts==points)));
        
        tasklabs = mainSig(vv).BoundLabels(inds2use,:);
        
        % cross validate leaving out one run at a time
%         cvlabs = runlabs;
        cvlabs = mainSig(vv).SessLabels(inds2use);
        nCV = numel(unique(cvlabs));
        
        predlabs = nan(size(dat,1),1);

        for tt = 1:nTasks
            for pp=1:nPairs

                pairinds = (ismember(pointlabs,points(pairs(pp,1),:),'rows') |...
                    ismember(pointlabs,points(pairs(pp,2),:),'rows')) & ...
                    tasklabs==tt;
                predlabs(pairinds) = nan;
                
                for cv=1:nCV

                    trninds = pairinds & cvlabs~=cv;
                    tstinds = pairinds & cvlabs==cv;

                    trnlabs = stim_inds(trninds);
                    assert(numel(unique(trnlabs))==2);
                    assert(sum(trnlabs==pairs(pp,1))==sum(trnlabs==pairs(pp,2)));

                    [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);

                    predlabs(tstinds) = label;

                end

                assert(~any(isnan(predlabs(pairinds))))

                acc = mean(predlabs(pairinds)==stim_inds(pairinds));
                fprintf('S%s, %s, task %d, pair %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},tt,pp, acc*100);
                all_acc(ss,vv,tt,pp) = acc;
                acc_matrix(ss,vv,tt,pairs(pp,1),pairs(pp,2)) = acc;
                acc_matrix(ss,vv,tt,pairs(pp,2),pairs(pp,1)) = acc;

            end
        end
    end
    
    tt=nTasks+1;
    
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    for vv=1:nVOIs

        inds2use = repSig(vv).IsMainLabels==1;

        dat = repSig(vv).dat_avg(inds2use,:);

        pointlabs = repSig(vv).PointLabels(inds2use,:);
        [pts, ~, stim_inds] = unique(pointlabs,'rows');
        assert(all(all(pts==points)));
      
        % cross validate leaving out one run at a time
        cvlabs = repSig(vv).SessLabels(inds2use);
        nCV = numel(unique(cvlabs));
        
        predlabs = nan(size(dat,1),1);

        
        for pp=1:nPairs

            pairinds = (ismember(pointlabs,points(pairs(pp,1),:),'rows') |...
                ismember(pointlabs,points(pairs(pp,2),:),'rows'));
            
            predlabs(pairinds) = nan;

            for cv=1:nCV

                trninds = pairinds & cvlabs~=cv;
                tstinds = pairinds & cvlabs==cv;

                trnlabs = stim_inds(trninds);
                assert(numel(unique(trnlabs))==2);
                assert(sum(trnlabs==pairs(pp,1))==sum(trnlabs==pairs(pp,2)));

                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),stim_inds(trninds));

                predlabs(tstinds) = label;

            end

            assert(~any(isnan(predlabs(pairinds))))

            acc = mean(predlabs(pairinds)==stim_inds(pairinds));
            fprintf('S%s, %s, repeat task, pair %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},pp, acc*100);
            all_acc(ss,vv,tt,pp) = acc;
            acc_matrix(ss,vv,tt,pairs(pp,1),pairs(pp,2)) = acc;
            acc_matrix(ss,vv,tt,pairs(pp,2),pairs(pp,1)) = acc;

        end
       
    end
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)

%% plot average pairwise decoding, one task per plot
for tt=1:nTasks+1
    mean_acc =  mean(mean(all_acc(:,:,tt,:),4),1);
    se_acc = squeeze(std(all_acc(:,:,tt,:),[],4))./sqrt(nPairs);

    plot_barsAndStars(mean_acc',se_acc',[],[],0.5,[0.4, 1],ROI_names,[],'Accuracy',sprintf('%s\nAverage of all pairwise classifiers',task_names{tt}),task_colors(tt,:))
end

%% plot average pairwise decoding, tasks side by side
v2plot = [1:5,10,11,6:9,12:17];

mean_acc = squeeze(mean(mean(all_acc(:,v2plot,:,:),4),1));

se_acc = squeeze(mean(std(all_acc(:,v2plot,:,:),[],4),1))./sqrt(nPairs);

plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names,'Accuracy','Average of all pairwise classifiers',task_colors)

set(gcf,'Position',[200,200,1000,400])

%% plot decoding versus distance in featurespace

v2plot = [1:5,10];
% v2plot =[12:17];
ss=1;

cols = plasma(numel(v2plot)+1);
cols = cols(1:numel(v2plot),:);

for tt=1:nTasks+1

    figure;hold all;

    for vv=1:numel(v2plot)


        meanvals = zeros(size(unvals));
        sevals = zeros(size(unvals));

        for uu=1:length(unvals)

            inds = pairwise_dist_round==unvals(uu);
            meanvals(uu) = mean(all_acc(ss,v2plot(vv),tt,inds));
            sevals(uu) = std(all_acc(ss,v2plot(vv),tt,inds))./sqrt(numel(inds));

        end

        errorbar(unvals, meanvals, sevals,'Color',cols(vv,:),'LineWidth',2);

        xlabel('Distance between shapes in feature space')
        ylabel('Decoding accuracy');


    end
    line(get(gca,'XLim'),[0.5, 0.5],'Color','k')

    title(sprintf('%s Task\nDecoding accuracy vs. distance in feature space',task_names{tt}));
    legend(ROI_names(v2plot),'Location','EastOutside');
    set(gcf,'Color','w');

    ylim([0.2, 1])

end



%% plot decoding along each axis separately
% separating out classifiers that go across the boundary versus same side.

v2plot = [1:5,10,11,6:9,12:17];

axis_names = {'Decode along dim 1','Decode along dim 2'};

across_names = {'across boundary', 'same side of boundary'};
across_vals = [1, 0];

for dd=1:2
    
    for aa=1:2
    
        pairs2use = which_dim==dd & across_bound(:,dd)==across_vals(aa);
        sum(pairs2use)
        
        mean_vals = zeros(nSubj, numel(v2plot), 2);
        se_vals = zeros(nSubj, numel(v2plot), 2);

        for tt=1:2
            mean_vals(:,:,tt) = squeeze(mean(all_acc(:,v2plot,tt,pairs2use),4));
            se_vals(:,:,tt) = squeeze(std(all_acc(:,v2plot,tt,pairs2use),[],4))./sqrt(sum(pairs2use));        
        end

    %     if nSubj==1
        mean_acc = squeeze(mean_vals);
        se_acc = squeeze(se_vals);
    %     else
    %         mean_acc =  mean(mean_vals,1);
    %         se_acc = std(mean_vals,[],1)./sqrt(nSubj);
    %     end
        plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names(1:2),...
            'Accuracy',sprintf('%s: %s',axis_names{dd}, across_names{aa}),task_colors(1:2,:));
        set(gcf,'Position',[200,200,1000,400])
    end
end   

%% visualizing what the above decoding is supposed to do...
aa=2;
dd=1;
figure;hold all;
pairs2use=find(which_dim==dd & across_bound(:,dd)==across_vals(aa));

for pp=1:numel(pairs2use)
    pt1 = points(pairs(pairs2use(pp),1),:);
    pt2 = points(pairs(pairs2use(pp),2),:);
    plot([pt1(1), pt2(1)], [pt1(2), pt2(2)], '-o')
end
set(gca,'XLim',[0,5]);
set(gca,'YLim',[0,5]);
plot(center_pos, [0,5],'k')
plot([0,5],center_pos,'k')
title(sprintf('%s: %s',axis_names{dd}, across_names{aa}))
%% plot decoding along each axis separately -overlay tasks
v2plot = [1:5,10,11,6:9,12:17];
axis_names = {'Decode along dim 1','Decode along dim 2'};
for dd=1:2
    
    pairs2use = which_dim==dd;
    sum(pairs2use)
    mean_vals = zeros(nSubj, numel(v2plot), 2);
    se_vals = zeros(nSubj, numel(v2plot), 2);
    
    for tt=1:2
        mean_vals(:,:,tt) = squeeze(mean(all_acc(:,v2plot,tt,pairs2use),4));
        se_vals(:,:,tt) = squeeze(std(all_acc(:,v2plot,tt,pairs2use),[],4))./sqrt(sum(pairs2use));        
    end

%     if nSubj==1
    mean_acc = squeeze(mean_vals);
    se_acc = squeeze(se_vals);
%     else
%         mean_acc =  mean(mean_vals,1);
%         se_acc = std(mean_vals,[],1)./sqrt(nSubj);
%     end
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names(1:2),'Accuracy',sprintf('%s',axis_names{dd}),task_colors(1:2,:));
    set(gcf,'Position',[200,200,1000,400])
end   

%% plot confusion matrix...sorting based on quadrants...
close all
v2plot = [2,10,5:9];

[~,sort_order] = sort(quadrant);

for vv=1:numel(v2plot)

    figure;hold all;
        
    for tt=1:nTasks+1

        subplot(2,2,tt);hold all;
        
        mat = squeeze(mean(acc_matrix(:,v2plot(vv),tt,:,:),1));
        mat_sorted = mat(sort_order,:);
        mat_sorted = mat_sorted(:,sort_order);
        
        imagesc(mat_sorted)
        
        axis('square')
        xlim([0.5, nStimsMain+0.5]) 
        ylim([0.5, nStimsMain+0.5]) 
        xticks([2.5:4:16]); xticklabels({'Q1','Q2','Q3','Q4'})
        yticks([2.5:4:16]); yticklabels({'Q1','Q2','Q3','Q4'})
        title(sprintf('%s',task_names{tt}));
        set(gca,'clim',[0,1])
        
    end
    
    suptitle(sprintf('%s',ROI_names{v2plot(vv)}));
    set(gcf,'Position',[200,200,1000,1000])
    set(gcf,'Color','w')
end


%% plot confusion matrix...sorting based on x-coordinate...
close all
% v2plot = [1:5,10];
v2plot = [2,10,5];

[~,sort_order] = sort(quadrant);
ptlabs = [];
for pp=1:size(points,1)
    ptlabs{pp} = sprintf('%.1f, %.1f',points(pp,1),points(pp,2));
end

for vv=1:numel(v2plot)

    figure;hold all;
        
    for tt=1:nTasks+1

        subplot(2,2,tt);hold all;
        
        mat = squeeze(mean(acc_matrix(:,v2plot(vv),tt,:,:),1));
       
        imagesc(mat)
        
%         mat_sorted = mat(sort_order,:);
%         mat_sorted = mat_sorted(:,sort_order);
        
%         imagesc(mat_sorted)
        
        axis('square')
        xlim([0.5, nStimsMain+0.5]) 
        ylim([0.5, nStimsMain+0.5]) 
        xticks([]);
        yticks([1:16]); yticklabels(ptlabs)
        title(sprintf('%s',task_names{tt}));
        set(gca,'clim',[0,1])
        
    end
    
    suptitle(sprintf('%s',ROI_names{v2plot(vv)}));
    set(gcf,'Position',[200,200,1000,1000])
    set(gcf,'Color','w')
end

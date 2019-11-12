%% analyze behavior of all subjects during 3 tasks

clear
close all

root = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/BehavioralPilot/';
addpath('/usr/local/serenceslab/maggie/mFiles/Classifiers/');

sublist = {'01','02'};
nSubj = length(sublist);

tasklist = {'ShapeBinaryTask1','ShapeBinaryTask2','ShapeQuadTask'};
tasklistplot = {'Binary (1)','Binary (2)','Quadrant'};
nTasks = length(tasklist);

nRunsEach = 5;

%% get full grid for binary tasks

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  
[gridx,gridy] = meshgrid(all_pts,all_pts);
all_grid_points = [gridx(:),gridy(:)];

% create psychometric curves - how likely is it that they report to the
% right/left of center boundary? [subj,task,position along relevant axis]
nRespRight = zeros(nSubj,2,numel(all_pts),nRunsEach);
nRespLeft = zeros(nSubj,2,numel(all_pts),nRunsEach);
nRespFail = zeros(nSubj,2,numel(all_pts),nRunsEach);

nDiffLevels = 13;

acclist = [];
difflist = [];
slist = [];
rlist = [];
tlist = [];

nrespgrid = zeros(nSubj, nTasks, numel(all_pts),numel(all_pts), 5);

%% load data

for ss=1:nSubj
    start_1 = [];stop_3 = [];
    for tt = 1:nTasks
       
        file = dir(fullfile(root, 'Data_V1',sprintf('*S%s*%s*.mat',sublist{ss},tasklist{tt})));
        if tt==1
            train_file = file(contains({file.name},'TRAINING'));
            load(fullfile(train_file.folder, train_file.name));
            start_1 = TheData(1).t.StartTime;
        end
        file = file(~contains({file.name},'TRAINING'));
        
        load(fullfile(file.folder, file.name));
        
        fprintf('Subject S%s did %d runs of %s\n',sublist{ss},length(TheData),tasklist{tt});
        
       
        if tt==3
            stop_3 = TheData(end).t.EndTime;
        end
       
        nTrials = numel(TheData(1).data.Response);
        
        for rr=1:nRunsEach
            
            acclist = [acclist; TheData(rr).data.Accuracy];
            difflist = [difflist; TheData(rr).p.difficulty];
            slist = [slist; ss];
            rlist = [rlist; rr];
            tlist = [tlist; tt];
            
%             dif = TheData(rr).p.difficulty;
%             accByDiff(ss,tt,dif,rr) = TheData(rr).data.Accuracy;
%             
            if tt<3
                ax = TheData(rr).p.which_bound;
                for tr = 1:nTrials
                    ind = find(all_pts==TheData(rr).p.imcoords(tr,ax));
                    if TheData(rr).data.Response(tr)==2
                        nRespRight(ss,tt,ind,rr) = nRespRight(ss,tt,ind,rr) + 1;
                    elseif TheData(rr).data.Response(tr)==1
                        nRespLeft(ss,tt,ind,rr) = nRespLeft(ss,tt,ind,rr) + 1;
                    else
                        nRespFail(ss,tt,ind,rr) = nRespFail(ss,tt,ind,rr) + 1;
                    end
                end
                assert(sum(nRespFail(ss,tt,:,rr))+sum(nRespRight(ss,tt,:,rr))+sum(nRespLeft(ss,tt,:,rr))==nTrials)
            end
            
            for tr =1:nTrials
                ind1 = find(all_pts==TheData(rr).p.imcoords(tr,1));
                ind2 = find(all_pts==TheData(rr).p.imcoords(tr,2));
                respind = TheData(rr).data.Response(tr);
                if respind==0
                    respind=5;
                end
                nrespgrid(ss,tt,ind1, ind2,respind) = nrespgrid(ss,tt,ind1, ind2,respind) + 1;
                
            end
            
        end
        
    end
    
    fprintf('subject %d did the expt in around %.2f minutes\n',ss,(stop_3-start_1)/60);
    
end

%% plot accuracy at each task
figure();hold all;

vals = zeros(nTasks, nSubj);

for ss=1:nSubj
       
    for tt=1:3
        
        yvals = acclist(slist==ss & tlist==tt);
        
        vals(tt,ss) = mean(yvals);
        
        
    end
    plot(1:3, vals(:,ss),'o');

end

meanvals = squeeze(nanmean(vals,2));
sevals = squeeze(nanstd(vals,[],2));
errorbar(1:3, meanvals,sevals,'Color','k');

xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Accuracy, averaged over runs'))

set(gca,'XLim',[0,4], 'YLim',[0.0, 1]);
set(gca,'XTick',1:3,'XTickLabel',tasklistplot,'XTickLabelRotation',45)

set(gcf,'Color','w');


%% plot accuracy versus difficulty level

for tt=1:3

    figure();hold all;

    for ss=1:nSubj
        
        xvals = difflist(slist==ss & tlist==tt);
        yvals = acclist(slist==ss & tlist==tt);
        
        plot(xvals, yvals, '-o')
        
        xlabel(sprintf('Difficulty'))
        ylabel(sprintf('Accuracy'));
        title(sprintf('Accuracy versus difficulty\n %s task',tasklistplot{tt}))
        
        set(gca,'XLim',[5,13], 'YLim',[0.0, 1]);
        
        if tt<3
            plot(get(gca,'XLim'), [0.5,0.5],'k')
        else
            plot(get(gca,'XLim'), [0.25,0.25],'k')
        end
    end
    set(gcf,'Color','w');
end

%%
allax=  [];
tick_pts = [1,24,47];

for tt=1:3
    
    figure;hold all;
    
    sumresp = squeeze(sum(nrespgrid(:,tt,:,:,:),1));
    probrespeach = sumresp(:,:,1:4)./sum(sumresp(:,:,1:4),3);
    never_shown = sum(sumresp(:,:,1:4),3)==0;
    never_shown = never_shown(:);
    
    npts = size(probrespeach,1);
    probrespeach = reshape(probrespeach, [size(probrespeach,1)*size(probrespeach,2), size(probrespeach,3)]);
    probrespeach(never_shown,:) = nan;
    
    probrespeach = reshape(probrespeach, [sqrt(size(probrespeach,1)),sqrt(size(probrespeach,1)), size(probrespeach,2)]);
    
    if tt<3
        for rr=1:2
            allax = [allax, subplot(1,2,rr)];
            sanePColor(1:npts,1:npts,probrespeach(:,:,rr));
            colorbar()
            title(sprintf('probability of resp %d',rr))
            axis square 
            set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts));
            set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YDir','normal');
        end
    else
        for rr=1:4
            allax = [allax, subplot(2,2,rr)];
            sanePColor(1:npts,1:npts,probrespeach(:,:,rr));
            colorbar()
            title(sprintf('probability of resp %d',rr))
            axis square 
            set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts));
            set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YDir','normal');
        end
    end
    
    suptitle(sprintf('%s task',tasklistplot{tt}))
    
    set(gcf,'Color','w')
end

match_clim(allax)
%% plot psychometric curves (not really enough data for this to look good)

 for tt=1:2

    figure();hold all;

    props = squeeze(nRespRight(:,tt,:,:)./(nRespLeft(:,tt,:,:)+nRespRight(:,tt,:,:)));
    props = squeeze(nanmean(props,3));

    meanvals = squeeze(nanmean(props,1));
    sevals = squeeze(nanstd(props,[],1));
    errorbar(all_pts, meanvals,sevals);


    xlabel(sprintf('Position along axis %d',tt))
    ylabel(sprintf('Proportion responded right'));
    title(sprintf('Psychometric curve, mean over runs\n %s task',tasklist{tt}))
%     xlim([-1.3,1.3])
% 
    plot([2.5,2.5], get(gca,'YLim'),'k')
% 
    plot(get(gca,'XLim'), [0.5,0.5],'k')

end
%% get full grid for binary tasks

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  
[gridx,gridy] = meshgrid(all_pts,all_pts);
all_grid_points = [gridx(:),gridy(:)];

which_bound = 1;
% now taking out images at exactly the prototype locations so that we
% never use these during task
nsteps_main = 4;
main_pts = round(linspace(start,stop, nsteps_main),1);
proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
if which_bound==1
    proto_coords = [proto_pts(2),center; proto_pts(1),center];
else
    proto_coords = [center, proto_pts(2); center, proto_pts(1)];
end
%     proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
assert(numel(proto_inds)==2)
all_grid_points(proto_inds,:) = [];

% Also taking out any images along boundary because these
% are ambiguous 
bound_inds = find(all_grid_points(:,which_bound)==center);
all_grid_points(bound_inds,:) = [];

% Next, for each point in the full grid, define the difficulty level 
% based on distance from the boundary
dist_from_bound = round(abs(all_grid_points(:,which_bound)-center),1);

% bin these values into 13 "levels"
[undist, ~ , dist_groups] = unique(round(dist_from_bound,1)); 

% define the start of each bin
bin_edges = round(fliplr([0.1:0.1:0.9,1.2:0.3:2.1]),1); 
dist_groups_binned = zeros(size(dist_groups));
mean_dist_bin_binary = zeros(numel(bin_edges),1);
for bb=1:numel(bin_edges)
    if bb>1
        inds = dist_from_bound>=bin_edges(bb) & dist_from_bound<bin_edges(bb-1);
    else
        inds = dist_from_bound>=bin_edges(bb);
    end
    dist_groups_binned(inds) = bb;
    mean_dist_bin_binary(bb) = mean(dist_from_bound(inds));
    
end  

nperbin_binary = sum(repmat(dist_groups_binned, 1, size(bin_edges,2))==repmat(1:size(bin_edges,2),size(dist_groups_binned,1),1),1);

%% get full grid for quadrant task

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  
[gridx,gridy] = meshgrid(all_pts,all_pts);
all_grid_points = [gridx(:),gridy(:)];

% now taking out images at exactly the prototype locations so that we
% never use these during task
nsteps_main = 4;
main_pts = round(linspace(start,stop, nsteps_main),1);
proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
assert(numel(proto_inds)==4)
all_grid_points(proto_inds,:) = [];

% Also taking out any images along quadrant boundaries because these
% are ambiguous 
bound_inds = find(any(all_grid_points==center,2));
all_grid_points(bound_inds,:) = [];

% Next, for each point in the full grid, define the difficulty level 
% based on distance from the nearest quadrant boundary
dist_from_bound = round(min(abs((all_grid_points-repmat(center, size(all_grid_points,1),2))),[],2),1);

% bin these values into 13 "levels"
[undist, ~ , dist_groups] = unique(round(dist_from_bound,1)); 

% define the start of each bin
bin_edges = round(fliplr([0.1:0.1:0.9,1.2:0.3:2.1]),1); 
dist_groups_binned = zeros(size(dist_groups));
mean_dist_bin_quadrant = zeros(numel(bin_edges),1);

for bb=1:numel(bin_edges)
    if bb>1
        inds = dist_from_bound>=bin_edges(bb) & dist_from_bound<bin_edges(bb-1);
    else
        inds = dist_from_bound>=bin_edges(bb);
    end
    dist_groups_binned(inds) = bb;
    mean_dist_bin_quadrant(bb) = mean(dist_from_bound(inds));
end  

nperbin_quadrant = sum(repmat(dist_groups_binned, 1, size(bin_edges,2))==repmat(1:size(bin_edges,2),size(dist_groups_binned,1),1),1);

%%

diff_to_dist = [mean_dist_bin_binary, mean_dist_bin_binary, mean_dist_bin_quadrant];

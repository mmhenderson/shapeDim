%% analyze behavior of all subjects during 3 tasks

clear
close all

root = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/BehavioralPilot/';
addpath('/usr/local/serenceslab/maggie/mFiles/Classifiers/');

% first list ALL subjects (even partial/incomplete), along with which 
% counter-balance condition they were in
sublist = [1:9];
cb_order = [1,2,1,2,1,1,2,2,1];

% now list any that I want to remove
sub2remove = [];
inds2remove = ismember(sublist,sub2remove);
sublist(inds2remove) = [];
cb_order(inds2remove) = [];

% how many are left?
nSubj = length(sublist);

% these are the strings that will be in the data files names for each task
tasklist = {'ShapeBinaryTaskV2_bound1','ShapeBinaryTaskV2_bound2','ShapeQuadTaskV2'};
taskorder1 = [1,1,2,2,3,3]';
taskorder2 = [2,2,1,1,3,3]';
maporder1 = [1,2,1,2,1,2]';
maporder2 = [1,2,1,2,3,4]';
% these are nicer strings for plotting
tasklistplot = {'Binary (1)','Binary (2)','Quadrant'};
nTasks = length(tasklist);

allacc = nan(nSubj, nTasks, 8);
acc_by_time = nan(nSubj, 2, 6, 2);
legendlabs = [];

plotRTs = 0;
plotAccByTime = 0;
plot2Dchoice=1;

%% create point grid in feature space

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  
[gridx,gridy] = meshgrid(all_pts,all_pts);
all_grid_points = [gridx(:),gridy(:)];

distances_from_bound = [-0.8:0.1:-0.1, 0.1:0.1:0.8];
nDistances = numel(distances_from_bound);

% create psychometric curves - how likely is it that they report to the
% right/left of center boundary? [subj,task,position along relevant axis]
nRespRight = zeros(nSubj,2,nDistances);
nRespLeft = zeros(nSubj,2,nDistances);
nRespFail = zeros(nSubj,2,nDistances);

nRespGrid = zeros(nSubj, nTasks, numel(all_pts),numel(all_pts), 5);

%% load data

for si = 1:nSubj

    subinit = sprintf('S%02.f', sublist(si));
    
    % look at all data in the folder, figure out the dates of the two
    % sessions
    alldat = dir(fullfile(root, 'Data_V2',sprintf('%s*.mat',subinit)));
    names = {alldat.name};
    dates = cellfun(@(x)x(end-9:end-4),names,'UniformOutput',false);
    sess_dates = unique(dates);
    [~,order] = sort(str2double(sess_dates),'ascend');
    sess_dates = sess_dates(order);
    nSess = numel(sess_dates);
    if nSess==1
        fprintf('%s: only one session\n',subinit)
    end
    % figure out which counter-balance condition this subject was in
    if cb_order(si)==1
        taskorders = [taskorder1, taskorder2];
        maporders = [maporder1, maporder2];
    else
        taskorders = [taskorder2, taskorder1];
        maporders = [maporder2, maporder1];
    end
    
    legendlabs{si} = sprintf('%s - cb %d',subinit,cb_order(si));
    
    task_counts = zeros(nTasks,1);
    
    allrts = [];
    zz = 0;
    for se = 1:nSess
        
        last_time = 0;
        
        % looping through tasks in order of when they were done in the
        % session. Note the tasks go in different orders for differnet
        % sessions
        for xx = 1:size(taskorders,1)
            
            curr_task = taskorders(xx,se);
            curr_map = maporders(xx,se);
           
            file = dir(fullfile(root, 'Data_V2',sprintf('%s_%s_mapping%d_%s.mat',subinit,tasklist{curr_task},curr_map,sess_dates{se})));
            
            if numel(file)==0
                file = dir(fullfile(root, 'Data_V2',sprintf('%s_%s_mapping%d_%s.mat',subinit,tasklist{curr_task},curr_map,sess_dates{3-se})));
                if numel(file)==0
                    fprintf('missing file for %s: task %d mapping %d session %d\n',subinit, curr_task, curr_map, se);
                    continue
                else
                    fprintf('%s:task %d mapping %d session %d was done in session %d instead\n',subinit, curr_task, curr_map, se,3-se);
                end
            end
   
            load(fullfile(file.folder, file.name));
            
            % make sure the runs were all done in the expected order, so
            % this one started after the previous one in the list
            if ~(str2double(TheData(1).t.TimeStamp)>last_time)
                fprintf('%s session %d: task %d mapping %d happened early\n',subinit, se,curr_task,curr_map);
            end
            last_time = str2double(TheData(1).t.TimeStamp);
   
            nTrials = numel(TheData(1).data.Response);
            
            nRuns = numel(TheData);
            
            for rr=1:nRuns
                
                % just a sanity check to make sure this file contains the
                % correct task
                if curr_task<3
                    assert(TheData(rr).p.which_bound==curr_task)
                else
                    assert(numel(TheData(rr).p.proto_order)==4)
                end
                
                % keep track of how many runs of the current task i have 
                task_counts(curr_task) = task_counts(curr_task) + 1;
                
                timeout = TheData(rr).data.Response==0 & ismember(TheData(rr).data.Response,TheData(rr).p.proto_order);
                
                if sum(timeout)>5
                    % print a warning if there were lots of timeouts
                    fprintf('%s sess %d task %d map %d run %d: %d timeout trials\n',subinit,se,curr_task,curr_map,rr,sum(timeout));
                end
                
                rts = TheData(rr).t.RespTimeFromOnset;
                rts = rts(~timeout);
                allrts = [allrts; rts];
                
                resp = TheData(rr).data.Response(~timeout);
                correct_resp = TheData(rr).p.category(~timeout);
                acc = mean(resp==correct_resp);
                
                allacc(si,curr_task,task_counts(curr_task)) = acc;
                acc_by_time(si,se,xx,rr) = acc;
                
                
                %% now look at trialwise performance, to make psychometric curves and such
                
                if curr_task<3
                    % make psychometric curves for the binary tasks only
                    for tr = 1:nTrials
%                         zz=zz+1
                        % for this trial, how far was the stimulus from
                        % boundary (difficulty)?
                        sign = double(TheData(rr).p.imcoords(tr,curr_task)>center);
                        sign(sign==0) = -1;
                        ind = find(round(distances_from_bound,1)==sign*round(TheData(rr).p.difficulty(tr),1));
                        assert(numel(ind)==1)
                        % what was their response on this trial?
                        resp = TheData(rr).data.Response(tr);
                        if resp~=0 && ismember(resp,TheData(rr).p.proto_order)
                            % resp is which button they pressed,
                            % resp_unmapped is which category this
                            % corresponded to in the active response
                            % mapping
                            resp_unmapped = TheData(rr).p.proto_order(resp);
                        else
                            resp_unmapped=0;
                        end
                        
                        % add to a running counter
                        if resp_unmapped==2
                            nRespRight(si,curr_task,ind) = nRespRight(si,curr_task,ind) + 1;
                        elseif resp_unmapped==1
                            nRespLeft(si,curr_task,ind) = nRespLeft(si,curr_task,ind) + 1;
                        else
                            nRespFail(si,curr_task,ind) = nRespFail(si,curr_task,ind) + 1;
                        end
                        
%                         nansum(nRespFail(si,curr_task,:))+nansum(nRespRight(si,curr_task,:))+nansum(nRespLeft(si,curr_task,:))
                        
                    end
                    assert(nansum(nRespFail(si,curr_task,:))+nansum(nRespRight(si,curr_task,:))+nansum(nRespLeft(si,curr_task,:))==nTrials*task_counts(curr_task))

                end
                
                % calculating prob of each choice based on where we are in
                % feature space, for all tasks
                for tr =1:nTrials
                    ind1 = find(all_pts==TheData(rr).p.imcoords(tr,1));
                    ind2 = find(all_pts==TheData(rr).p.imcoords(tr,2));
                    
                    % what was their response on this trial?
                    resp = TheData(rr).data.Response(tr);
                    if resp~=0 && ismember(resp,TheData(rr).p.proto_order)
                        % resp is which button they pressed,
                        % resp_unmapped is which category this
                        % corresponded to in the active response
                        % mapping
                        resp_unmapped = TheData(rr).p.proto_order(resp);
                    else
                        resp_unmapped=0;
                    end
                    
                    respind = resp_unmapped;
                    if respind==0
                        respind=5;
                    end
                    nRespGrid(si,curr_task,ind1,ind2,respind) = nRespGrid(si,curr_task,ind1,ind2,respind) + 1;

                end
            end
        end
    end
    %%
    if plotRTs
    figure;hold all;
    scatter(1:numel(allrts),allrts)
    title(sprintf('all RTs for %s',subinit));
    ylim([0,2]);
    xlabel('trials in whole expt')
    end
    %% check a few things

    for tt=1:nTasks
        if tt<3
            assert(nansum(nRespFail(si,tt,:))+nansum(nRespRight(si,tt,:))+nansum(nRespLeft(si,tt,:))==nTrials*task_counts(tt))
        end
        assert(nansum(nansum(nansum(nRespGrid(si,tt,:,:,:))))==nTrials*task_counts(tt))
    end
end


%% plot accuracy at each task
figure();hold all;

vals = nanmean(allacc, 3);

for si=1:nSubj
   
    plot(1:3, vals(si,:),'-o');

end

meanvals = squeeze(nanmean(vals,1));
sevals = squeeze(nanstd(vals,[],1))./sqrt(nSubj);
errorbar(1:3, meanvals,sevals,'Color','k');
plot([0,2.5],[0.5, 0.5],'Color','k','LineStyle','--');
plot([2.5,4],[0.25, 0.25],'Color','k','LineStyle','--');
xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Accuracy, averaged over runs'))

set(gca,'XLim',[0,4], 'YLim',[0.0, 1]);
set(gca,'XTick',1:3,'XTickLabel',tasklistplot,'XTickLabelRotation',45)

set(gcf,'Color','w');

legend(legendlabs,'Location','EastOutside')


%% plot accuracy in order of when blocks were done
if plotAccByTime
for se = 1:2
    
    for cb = 1:2

        figure();hold all;

        sub2use = find(cb_order==cb);
        
        for si=sub2use

            dat = squeeze(acc_by_time(si,se,:,:))';
            dat = dat(:);
            plot(1:numel(dat), dat,'-o');

        end

        xlabel(sprintf('Task'))
        ylabel(sprintf('Accuracy'));
        title(sprintf('Session %d, cb %d\naccuracy versus run number',se,cb))

        ylim([0,1]);
        xlim([0, numel(dat)+1]);

        linepos = [4.5, 8.5];
        for ll=1:numel(linepos)
            line([linepos(ll),linepos(ll)],get(gca,'YLim'),'Color','k');
        end

        set(gcf,'Color','w');

        legend(legendlabs(sub2use),'Location','EastOutside')
    end
end
end
%% plot choices as a function of where in feature space we are
if plot2Dchoice
    
allax=  [];
tick_pts = [1,24,47];

for tt=1:3
    
    figure;hold all;
    
    sumresp = squeeze(sum(nRespGrid(:,tt,:,:,:),1));
    probrespeach = sumresp(:,:,1:4)./sum(sumresp(:,:,1:4),3);
    never_shown = sum(sumresp(:,:,1:4),3)==0;
    never_shown = never_shown(:);
    
    npts = size(probrespeach,1);
    probrespeach = reshape(probrespeach, [size(probrespeach,1)*size(probrespeach,2), size(probrespeach,3)]);
    probrespeach(never_shown,:) = nan;
    
    probrespeach = reshape(probrespeach, [sqrt(size(probrespeach,1)),sqrt(size(probrespeach,1)), size(probrespeach,2)]);
    
    if tt<3
        for rr=1:2
            allax = [allax, subplot(1,2,rr)];hold all;
            % taking a transpose here because I find it confusing otherwise
            sanePColor(1:npts,1:npts,probrespeach(:,:,rr)');
%             plot(1, 2.5, 'ro')
            colorbar()
            title(sprintf('probability of resp %d',rr))
            axis square 
            set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts));
            set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YDir','normal');
            xlabel('feature 1');
            ylabel('feature 2');
        end
        
    else
        for rr=1:4
            allax = [allax, subplot(2,2,rr)];
            sanePColor(1:npts,1:npts,probrespeach(:,:,rr)');
            colorbar()
            title(sprintf('probability of resp %d',rr))
            axis square 
            set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts));
            set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YDir','normal');
            xlabel('feature 1');
            ylabel('feature 2')
        end
    end
    
    suptitle(sprintf('%s task',tasklistplot{tt}))
    
    set(gcf,'Color','w')
end

match_clim(allax)

end
%% plot psychometric curves (not really enough data for this to look good)

 for tt=1:2

    figure();hold all;
    props = squeeze(nRespRight(:,tt,:)./(nRespLeft(:,tt,:)+nRespRight(:,tt,:)));
    
    for si=1:nSubj
        
        plot(distances_from_bound,squeeze(props(si,:)),'-','Color',[0.8, 0.8, 0.8])
        
    end

    
    meanvals = squeeze(nanmean(props,1));
    sevals = squeeze(nanstd(props,[],1))./sqrt(nSubj);
    errorbar(distances_from_bound, meanvals,sevals,'Color','k');


    xlabel(sprintf('distance from center'))
    ylabel(sprintf('Proportion responded right'));
    title(sprintf('Psychometric curve, mean over runs\n %s task',tasklistplot{tt}))
    xlim([-1,1])
    ylim([0,1])
    plot([0,0], get(gca,'YLim'),'k')
% 
    plot(get(gca,'XLim'), [0.5,0.5],'k')

 end

 set(gcf,'Color','w')

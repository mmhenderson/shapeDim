%% analyze behavior of a single subject from their training day
% main goal is to see which tasks are hardest for them so that we can start
% to adjust the difficulty as soon as the scan session begins to try and
% match performance across all taks.

%%
clear
close all

currdir = pwd;
ndirsup=1;
filesepinds = find(currdir==filesep);
root = currdir(1:filesepinds(end-ndirsup+1));
addpath('/usr/local/serenceslab/maggie/mFiles/Classifiers/');

% first list the subject
sublist = [1];

% now list any that I want to remove
sub2remove = [];
inds2remove = ismember(sublist,sub2remove);
sublist(inds2remove) = [];

% how many are left?
nSubj = length(sublist);

% list the orders in which things will happen across entire experiment
% [nParts x nSess]
% which boundary is active? this defines which categorization dim they're using
bound_list_all = [1,2,3; 2,3,1; 3,1,2; 1,2,3; 2,3,1; 3,1,2];
% which response mapping? e.g. which finger for which category?
map_list_all = [1,2,1; 1,2,2; 1,2,1; 2,1,2; 2,1,1; 2,1,2];

% names of the tasks
tasklist = {'Linear (1)','Linear (2)','Checker'};
nTasks = length(tasklist);

nRunsPerPart = 2;
nPartsPerTask = 2;  % 2 mappings per task
nPartsTotal  = nTasks*nPartsPerTask;
nSess = 1;
acc_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
diff_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_by_run = nan(nSubj, nSess, nPartsTotal, nRunsPerPart);

legendlabs = [];

plotRTs = 1;
plotAccByTime = 1;
plot2DChoice=1;
plot1DChoice=1;
% 
%% quality control settings...
% what percent of trials can be timeout before i toss the whole block?
% if they fail to respond to 25 percent, probably they were confused and/or
% distracted and/or unmotivated.
timeout_pct_cutoff = 0.25;

%% create point grid in feature space

% first I'm defining all possible images that we can use in this task.
start = 0;    % min value along each axis
stop = 5; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = round(start:step:stop, 1);  

distances_from_bound = [-2.4, -0.8:step:-step, step:step:0.8, 2.4];
nDistances = numel(distances_from_bound);

% create psychometric curves - how likely is it that they report to the
% right/left of center boundary? [subj,task,position along relevant axis]
% for linear tasks only
nRespRight = zeros(nSubj,2,nDistances);
nRespLeft = zeros(nSubj,2,nDistances);
nRespFail = zeros(nSubj,2,nDistances);

nRespGrid = zeros(nSubj, nTasks, numel(all_pts),numel(all_pts), 3);

possible_responses =[1,2];
%% load data

actual_cb = zeros(nSubj,1);
for si = 1:nSubj

    subinit = sprintf('S%02.f', sublist(si));
    datadir = fullfile(root, 'DataBehavior',subinit,'Training');
    % look at all data in the folder, figure out the dates of the two
    % sessions
    alldat = dir(fullfile(datadir, sprintf('%s_MainTaskMRI_feedback*.mat',subinit)));
    names = {alldat.name};
    dates = cellfun(@(x)x(end-9:end-4),names,'UniformOutput',false);
    sess_dates = unique(dates);
    [~,order] = sort(str2double(sess_dates),'ascend');
    sess_dates = sess_dates(order);
    nSess = numel(sess_dates);
    se=1;
    
    % figure out which counter-balance condition this subject was in
    cb_ind = mod(sublist(si),3);
    actual_cb(si) = cb_ind;
    fprintf('%s: cb=%d\n',subinit, cb_ind);
    
    task_list = bound_list_all(:,cb_ind);
    map_list = map_list_all(:,cb_ind);
    position_list = map_list;
    
    legendlabs{si} = sprintf('%s',subinit);
    
    % counting tasks - combining both parts of each task together 
    task_counts = zeros(nTasks,1);
    
    allrts = [];

    % looping over the 6 parts of this session
    for xx = 1:nPartsTotal

        filename = fullfile(datadir,sprintf('%s_MainTaskMRI_feedback_sess%d_part%d_%s.mat',subinit,se,xx,sess_dates{se}));
        load(filename);
        p = TheData(1).p;

%             if se==1 && xx==1
%                 if ~mod(p.SubNum,2)
%                     fprintf('%s: even\n',subinit);
%                 else
%                     fprintf('%s: odd\n',subinit);
%                 end
%             end
%             fprintf('    %s: task=%d, map=%d, position=%d\n',sess_dates{se},p.which_bound, p.which_mapping, p.which_pos);
        % check and make sure the parameters in the file are correct
        assert(p.which_bound==task_list(xx,se));
        assert(p.which_mapping==map_list(xx,se));
        assert(p.which_pos==position_list(xx,se));

        curr_task = task_list(xx,se);
        curr_map = map_list(xx,se);

        nTrials = numel(TheData(1).data.Response);

%             if curr_map==2
%                 continue
%             end

        for rr=1:nRunsPerPart
            if size(TheData,2)<rr
                fprintf('Subject %s is missing run %d for sess %d part %d\n',p.SubNumStr,rr,se,xx);
                continue
            end

            p = TheData(rr).p;
            data = TheData(rr).data;
            t = TheData(rr).t;

            timeout = data.Response==0 | ~ismember(data.Response,possible_responses);

            if sum(timeout)>nTrials*timeout_pct_cutoff
                % print a warning if there were lots of timeouts
                fprintf('%s sess %d part %d run %d: %d timeout trials\n',subinit,se,xx,rr,sum(timeout));
                fprintf('Not using this block!\n');
                continue
            end

            % keep track of how many runs of the current task i have 
            task_counts(curr_task) = task_counts(curr_task) + 1;

            rts = t.RespTimeFromOnset;
%                 rts = rts(~timeout);

            allrts = [allrts, rts];

%                 resp = p.category(~timeout)';
            resp = data.Response(~timeout);
            correct_resp = p.category(~timeout);
            acc = mean(resp==correct_resp);

            acc_by_task(si,curr_task,task_counts(curr_task)) = acc;
            acc_by_run(si,se,xx,rr) = acc;
            diff_by_task(si,curr_task, task_counts(curr_task)) = p.RunDifficulty;
%                 cat_groups = [p.proto_order(1:2);p.proto_order(3:4)];

            % define which category was which - for mapping 2, this
            % would have been reversed relative to mapping 1, so make
            % them the same again.
            if curr_map==2
                category_unmapped = 3-p.category;
            else
                category_unmapped = p.category;
            end

            %% now look at trialwise performance, to make psychometric curves and such

            if curr_task<3
                % make psychometric curves for the linear tasks only
                for tr = 1:nTrials

                    % for this trial, how far was the stimulus from
                    % boundary (difficulty)?
                    sign = double(p.points(tr,curr_task)>center);
                    sign(sign==0) = -1;
                    dist = round(p.dist_from_bound(tr),1);

                    ind = find(round(distances_from_bound,1)==sign*dist,1);
                    assert(numel(ind)==1)
                    % what was their response on this trial?
                    if ~timeout(tr)
%                             resp = p.category(tr);
                        resp = data.Response(tr);                        
                        % resp_cat is which category their response
                        % indicated, which is swapped on different
                        % mappings. Make them the same again.
                        if curr_map==1
                            resp_unmapped = resp;
                        else
%                                 error('map should be 1')
                            resp_unmapped = 3-resp;
                        end
                    else
                        resp_unmapped=0;
                    end

                    % add to a running counter
                    if resp_unmapped==1
                        nRespRight(si,curr_task,ind) = nRespRight(si,curr_task,ind) + 1;
                    elseif resp_unmapped==2
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
                ind1 = find(all_pts==p.points(tr,1));
                ind2 = find(all_pts==p.points(tr,2));

                % what was their response on this trial?
                if ~timeout(tr)
%                         resp = p.category(tr);
                    resp = data.Response(tr);
                    % resp_cat is which category their response
                    % indicated, which is swapped on different
                    % mappings (the variable p.category has a different
                    % meaning).
                    if curr_map==1
                        resp_unmapped = resp;
                    else
                        resp_unmapped = 3-resp;
                    end
                else
                    resp_unmapped=0;
                end

                respind = resp_unmapped;
                if respind==0
                    respind=3;
                end
                nRespGrid(si,curr_task,ind1,ind2,respind) = nRespGrid(si,curr_task,ind1,ind2,respind) + 1;

            end
        end


        
    end
    %%
    run_colors = viridis(size(allrts,2));
    if plotRTs
        figure;hold all;
        ylim([0,2]);
        for rt = 1:size(allrts,2)
            xpts = (1:nTrials)+(rt-1)*nTrials;
            plot(xpts,allrts(:,rt),'o','Color',run_colors(rt,:));
            plot([xpts(end)+0.5,xpts(end)+0.5],get(gca,'YLim'),'-','Color',[0.8, 0.8, 0.8]); 
        end
        title(sprintf('all RTs for %s',subinit));
        set(gca,'XTick',nTrials/2:nTrials:nTrials*rt, 'XTickLabels',1:rt)
        xlabel('runs in whole expt')
        ylabel('RT (sec)')
        set(gcf,'Color','w')
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
sub_colors = viridis(nSubj+2);
vals = nanmean(acc_by_task, 3);

for si=1:nSubj
   
    plot(1:3, vals(si,:),'-o','Color',sub_colors(si,:));

end

meanvals = squeeze(nanmean(vals,1));
sevals = squeeze(nanstd(vals,[],1))./sqrt(nSubj);
errorbar(1:3, meanvals,sevals,'Color','k');
plot([0,4],[0.5, 0.5],'Color','k','LineStyle','--');
xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Accuracy, averaged over all runs'))

set(gca,'XLim',[0,4], 'YLim',[0.0, 1]);
set(gca,'XTick',1:3,'XTickLabel',tasklist,'XTickLabelRotation',45)

set(gcf,'Color','w');

legend(legendlabs,'Location','EastOutside')

%% print accuracy each task
vals = nanmean(acc_by_task, 3);
meanvals = squeeze(nanmean(vals,1));
meandiff = squeeze(nanmean(nanmean(diff_by_task,1),3));
fprintf('\nOverall performance:\n')
for tt = 1:nTasks
   fprintf('%s accuracy = %.2f pct, average diff = %.2f\n',tasklist{tt}, meanvals(tt)*100, meandiff(tt)) 
end
%% plot accuracy in order of when blocks were done
if plotAccByTime
    sub_colors = viridis(nSubj+2);
    chance_levels = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
for se = 1:nSess
  
    figure();hold all;

    sub2use = sublist;

    ylim([0,1]);
    xlim([0, nPartsTotal*nRunsPerPart+1]);
    lax = [];
    for si=sub2use

        dat = squeeze(acc_by_run(si,se,:,:));
        
        for pp = 1:nPartsTotal
            inds = nRunsPerPart*(pp-1)+1:nRunsPerPart*pp;
            myax = plot(inds, dat(pp,:),'-o','Color',sub_colors(si,:));
            if pp==1
                lax = [lax, myax];
            end
            plot(inds,chance_levels(pp)*ones(size(inds)),'--','Color',[0.8, 0.8, 0.8]);
        end
        
    end

    xlabel(sprintf('Run Number'))
    ylabel(sprintf('Accuracy'));
    title(sprintf('Session %d, all subjects\naccuracy versus run number',se))
    set(gca,'XTick',1:nPartsTotal*nTasks, 'XTickLabels',repmat((1:nPartsTotal)',nTasks,1))
   
    linepos = [nRunsPerPart+.5, nRunsPerPart*2+.5];
    for ll=1:numel(linepos)
        line([linepos(ll),linepos(ll)],get(gca,'YLim'),'Color','k');
    end

    set(gcf,'Color','w');

    legend(lax,legendlabs(sub2use),'Location','EastOutside')

end
end
%% plot choices as a function of where in feature space we are
if plot2DChoice
    
allax=  [];
tick_pts = [1,25,50];

for tt=1:3
    
    figure;hold all;
    
    sumresp = squeeze(sum(nRespGrid(:,tt,:,:,:),1));
    probrespeach = sumresp(:,:,1:2)./sum(sumresp(:,:,1:2),3);
    never_shown = sum(sumresp(:,:,1:2),3)==0;
    never_shown = never_shown(:);
    
    npts = size(probrespeach,1);
    probrespeach = reshape(probrespeach, [size(probrespeach,1)*size(probrespeach,2), size(probrespeach,3)]);
    probrespeach(never_shown,:) = nan;
    
    probrespeach = reshape(probrespeach, [sqrt(size(probrespeach,1)),sqrt(size(probrespeach,1)), size(probrespeach,2)]);

    for rr=1:2
        allax = [allax, subplot(1,2,rr)];hold all;
        % taking a transpose here because I find it confusing otherwise
        sanePColor(1:npts,1:npts,probrespeach(:,:,rr)');
        colormap(viridis)
%             plot(1, 2.5, 'ro')
        colorbar()
        title(sprintf('probability of resp %d',rr))
        axis square 
        set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts),'XLim',[tick_pts(1)-1, tick_pts(3)+1]);
        set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YLim',[tick_pts(1)-1, tick_pts(3)+1], 'YDir','normal');
        xlabel('feature 1');
        ylabel('feature 2');
        plot([tick_pts(2),tick_pts(2)], get(gca,'YLim'),'k');
        plot(get(gca,'XLim'),[tick_pts(2),tick_pts(2)], 'k');
    end
    
    suptitle(sprintf('%s task',tasklist{tt}))
    
    set(gcf,'Color','w')
end

match_clim(allax)

end
%% plot psychometric curves (only for the linear tasks)

cols = viridis(nSubj+2);

if plot1DChoice
    
    for tt=1:2
        
        figure();hold all;
        props = squeeze(nRespRight(:,tt,:)./(nRespLeft(:,tt,:)+nRespRight(:,tt,:)));
        
        for si=1:nSubj
            
            vals = squeeze(props(si,:));
            not_sampled = isnan(vals);
%             plot(distances_from_bound(~not_sampled),vals(~not_sampled),'-','Color',[0.8, 0.8, 0.8])
            plot(distances_from_bound(~not_sampled),vals(~not_sampled),'-','Color',cols(si,:))
           
            
        end
        
        if nSubj>1
            meanvals = squeeze(nanmean(props,1));
            sevals = squeeze(nanstd(props,[],1))./sqrt(sum(~isnan(props),1));
            not_sampled = isnan(meanvals);
            sevals = sevals(not_sampled);
        else
            meanvals = props;
            not_sampled = isnan(meanvals);
            sevals = [];
        end
        
        not_sampled = isnan(meanvals);
        
        errorbar(distances_from_bound(~not_sampled), meanvals(~not_sampled),sevals,'Color','k');
        
        
        xlabel(sprintf('distance from center'))
        ylabel(sprintf('Proportion responded right'));
        title(sprintf('Psychometric curve, mean over runs\n %s task',tasklist{tt}))
        xlim([distances_from_bound(1)-0.1, distances_from_bound(end)+0.1])
        ylim([0,1])
        plot([0,0], get(gca,'YLim'),'k')
        %
        plot(get(gca,'XLim'), [0.5,0.5],'k')
        set(gcf,'Color','w')
        
        legend(legendlabs,'Location','EastOutside')

    end
    
    
end
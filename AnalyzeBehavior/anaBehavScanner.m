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
sesslist = {[1,2,3]};

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

% info about main categorization tasks, initialized here
nRunsPerPart = 2;
nPartsPerTask = 2;  % 2 mappings per task
nPartsTotal  = nTasks*nPartsPerTask;
nSess = 3;
acc_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_easier_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_harder_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
diff_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_by_run = nan(nSubj, nSess, nPartsTotal, nRunsPerPart);
nTrialsPerRun=48;
allrts = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess, nTrialsPerRun);

% info about repeat detection task, initialize
nRunsRepeat = 4;
acc_repeat_task = nan(nSubj, nSess*nRunsRepeat);
acc_easier_repeat_task = nan(nSubj, nSess*nRunsRepeat);
acc_harder_repeat_task = nan(nSubj, nSess*nRunsRepeat);
diff_repeat_task = nan(nSubj, nSess*nRunsRepeat);
legendlabs = [];

%% choices for what to plot

plotRTs = 1;
plotAccByTime = 0;
plot2DChoice=0;
plot1DChoice=0;
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
    sess2do = sesslist{si};
    
    % counting tasks - combining both parts of each task together 
    task_counts = zeros(nTasks,1);
    repeat_task_count = 0;
    
    
    for se = 1:length(sess2do)
        
        sessstr = sprintf('Session%d',sess2do(se));
        datadir = fullfile(root, 'DataBehavior',subinit,sessstr);
        
        % look at all data in the folder, figure out the date of the sess
        alldat = dir(fullfile(datadir, sprintf('%s_MainTaskMRI_scannerversion*.mat',subinit)));
        sess_date = alldat(1).name(end-9:end-4);

        % figure out which counter-balance condition this subject was in
        cb_ind = mod(sublist(si),3);
        actual_cb(si) = cb_ind;
        fprintf('%s: cb=%d\n',subinit, cb_ind);

        task_list = bound_list_all(:,cb_ind);
        map_list = map_list_all(:,cb_ind);
        position_list = map_list;

        legendlabs{si} = sprintf('%s',subinit);

        
        allrts_this_sess = [];

        % looping over the 6 parts of this session
        for xx = 1:nPartsTotal

            filename = fullfile(datadir,sprintf('%s_MainTaskMRI_scannerversion_sess%d_part%d_%s.mat',subinit,sess2do(se),xx,sess_date));
            load(filename);
            p = TheData(1).p;

            % check and make sure the parameters in the file are correct
            assert(p.which_bound==bound_list_all(xx,sess2do(se)));
            assert(p.which_mapping==map_list_all(xx,sess2do(se)));
           
            curr_task = bound_list_all(xx,sess2do(se));
            curr_map = map_list_all(xx,sess2do(se));

            nTrials = numel(TheData(1).data.Response);

            for rr=1:nRunsPerPart
                if size(TheData,2)<rr
                    fprintf('Subject %s is missing run %d for sess %d part %d\n',p.SubNumStr,rr,sess2do(se),xx);
                    continue
                end

                p = TheData(rr).p;
                data = TheData(rr).data;
                t = TheData(rr).t;

                timeout = data.Response==0 | ~ismember(data.Response,possible_responses);

                if sum(timeout)>nTrials*timeout_pct_cutoff
                    % print a warning if there were lots of timeouts
                    fprintf('%s sess %d part %d run %d: %d timeout trials\n',subinit,sess2do(se),xx,rr,sum(timeout));
                    fprintf('Not using this block!\n');
                    continue
                end

                % keep track of how many runs of the current task i have 
                task_counts(curr_task) = task_counts(curr_task) + 1;

                rts = t.RespTimeFromOnset;
    %                 rts = rts(~timeout);

                allrts_this_sess = [allrts_this_sess, rts];

                allrts(si,curr_task,task_counts(curr_task),:) = rts;
                
    %                 resp = p.category(~timeout)';
                resp = data.Response(~timeout);
                correct_resp = p.category(~timeout);
                acc = mean(resp==correct_resp);

                acc_by_task(si,curr_task,task_counts(curr_task)) = acc;
                acc_easier_by_task(si,curr_task,task_counts(curr_task)) = data.MainGridAccuracy;
                acc_harder_by_task(si,curr_task,task_counts(curr_task)) = data.VariableAccuracy;
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
        run_colors = viridis(size(allrts_this_sess,2));
        if plotRTs
            figure;hold all;
            ylim([0,2]);
            for rt = 1:size(allrts_this_sess,2)
                xpts = (1:nTrials)+(rt-1)*nTrials;
                plot(xpts,allrts_this_sess(:,rt),'o','Color',run_colors(rt,:));
                plot([xpts(end)+0.5,xpts(end)+0.5],get(gca,'YLim'),'-','Color',[0.8, 0.8, 0.8]); 
            end
            title(sprintf('all RTs for %s, sess %d',subinit, sess2do(se)));
            set(gca,'XTick',nTrials/2:nTrials:nTrials*rt, 'XTickLabels',1:rt)
            xlabel('runs in whole sess')
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
        
        %% also loading data for repeat detection task (one-back)
        filename = fullfile(datadir,sprintf('%s_OneBackTaskMRI_sess%d_%s.mat',subinit,sess2do(se),sess_date));
        load(filename);
        for rr = 1:length(TheData)
            
            repeat_task_count = repeat_task_count+1;
            
            p = TheData(rr).p;
            data = TheData(rr).data;
            acc_repeat_task(si,repeat_task_count) = data.Accuracy;
            acc_easier_repeat_task(si,repeat_task_count) = data.MainGridAccuracy;
            acc_harder_repeat_task(si,repeat_task_count) = data.VariableAccuracy;
            diff_repeat_task(si,repeat_task_count) = p.RunDifficulty;
        end
    end
end


%% plot accuracy at each task
figure();hold all;
sub_colors = viridis(nSubj+2);
vals = nanmean(acc_by_task, 3);
vals = [vals, nanmean(acc_repeat_task,2)];
for si=1:nSubj
   
    plot(1:4, vals(si,:),'-o','Color',sub_colors(si,:));

end

meanvals = squeeze(nanmean(vals,1));
sevals = squeeze(nanstd(vals,[],1))./sqrt(nSubj);
errorbar(1:4, meanvals,sevals,'Color','k');
plot([0,4],[0.5, 0.5],'Color','k','LineStyle','--');
xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Accuracy, averaged over all runs'))

set(gca,'XLim',[0,5], 'YLim',[0.0, 1]);
set(gca,'XTick',1:4,'XTickLabel',[tasklist, 'One-Back'],'XTickLabelRotation',45)

set(gcf,'Color','w');

legend(legendlabs,'Location','EastOutside')

%% plot accuracy at each task - break down into easy and hard.

vals_easier = nanmean(acc_easier_by_task, 3);
vals_easier = [vals_easier, nanmean(acc_easier_repeat_task,2)];
vals_harder = nanmean(acc_harder_by_task, 3);
vals_harder = [vals_harder, nanmean(acc_harder_repeat_task,2)];
vals = nanmean(acc_by_task, 3);
vals = [vals, nanmean(acc_repeat_task,2)];

easyhard_colors = plasma(4);
easyhard_colors = easyhard_colors(2:3,:);
for si=1:nSubj
   
    figure();hold all;
    plot(1:4, vals_easier(si,:),'-o','Color',easyhard_colors(1,:));
    plot(1:4, vals_harder(si,:),'-o','Color',easyhard_colors(2,:));
    plot(1:4, vals(si,:),'-o','Color','k');
  
    plot([0,4],[0.5, 0.5],'Color','k','LineStyle','--');
    xlabel(sprintf('Task'))
    ylabel(sprintf('Accuracy'));
    subinit = sprintf('S%02.f', sublist(si));
    title(sprintf('Accuracy for %s\n averaged over all runs',subinit))

    set(gca,'XLim',[0,5], 'YLim',[0.0, 1]);
    set(gca,'XTick',1:4,'XTickLabel',[tasklist,'One-Back'],'XTickLabelRotation',45)

    set(gcf,'Color','w');

    legend({'easier trials','harder trials','all'},'Location','EastOutside')

end


%% print accuracy each task  

for si = 1:nSubj
    subinit = sprintf('S%02.f', sublist(si));
    meanvals = nanmean(acc_by_task(si,:,:), 3);  
    meanvals_easier = nanmean(acc_easier_by_task(si,:,:), 3);  
    meanvals_harder = nanmean(acc_harder_by_task(si,:,:), 3);  
    meandiff = nanmean(diff_by_task(si,:,:), 3);  
    fprintf('\n%s Overall performance:\n',subinit)
    for tt = 1:nTasks
       fprintf('    %s accuracy = %.2f pct, average diff = %.2f\n',tasklist{tt}, meanvals(tt)*100, meandiff(tt)) 
       fprintf('        easier = %.2f pct, harder = %.2f\n',meanvals_easier(tt)*100, meanvals_harder(tt)*100) 
    end
    fprintf('    Repeat task accuracy = %.2f pct, average diff = %.2f\n',nanmean(acc_repeat_task(si,:))*100, nanmean(diff_repeat_task(si,:)))
    fprintf('        easier = %.2f pct, harder = %.2f\n',nanmean(acc_easier_repeat_task(si,:))*100, nanmean(acc_harder_repeat_task(si,:))*100) 
end

meanvals = squeeze(nanmean(nanmean(acc_by_task, 3),1));
semvals = squeeze(nanstd(nanmean(acc_by_task,3),[],1))./sqrt(nSubj);
meandiff = squeeze(nanmean(nanmean(diff_by_task,3),1));
fprintf('\nOverall performance avg all subs:\n')
for tt = 1:nTasks
   fprintf('    %s accuracy = %.2f +/- %.2f pct, average diff = %.2f\n',tasklist{tt}, meanvals(tt)*100, semvals(tt)*100, meandiff(tt)) 
end
fprintf('    repeat task accuracy = %.2f +/- %.2f pct, average diff = %.2f\n',...
    nanmean(nanmean(acc_repeat_task,2),1)*100, nanstd(nanmean(acc_repeat_task,2),[],1)./sqrt(nSubj)*100,...
    nanmean(nanmean(diff_repeat_task,2),1))

%% plot accuracy in order of when blocks were done
if plotAccByTime
    sub_colors = viridis(nSubj+2);
    chance_levels = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    for se = 1:numel(sess2do)

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
                plot([inds(end)+0.5,inds(end)+0.5],get(gca,'YLim'),'-','Color',[0.8, 0.8, 0.8]); 
            end

        end

        xlabel(sprintf('Run Number'))
        ylabel(sprintf('Accuracy'));
        title(sprintf('Session %d, all subjects\naccuracy versus run number',sess2do(se)))
        set(gca,'XTick',1:nPartsTotal*nTasks, 'XTickLabels',repmat((1:nPartsTotal)',nTasks,1))
        
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

%% 
xlimhist_rt=[0, 3];
ylimhist_rt = [0,120];
task_colors = viridis(nTasks+1);
task_colors = task_colors(1:end-1,:);

if plotRTs
   
    figure;hold all;
    xx=0;
    for tt = 1:nTasks
       
            xx=xx+1;
            subplot(3,1,xx);hold all;
    %             figure;hold all;
            title(tasklist{tt})
            rt_vals = allrts(:,tt,:,:);
            rt_vals = rt_vals(:);
            fprintf('%d nans\n',sum(isnan(rt_vals)));
            rt_vals = rt_vals(~isnan(rt_vals));
            
            meanrt = mean(rt_vals);  
            stdev = std(rt_vals);
            assert(all(rt_vals<xlimhist_rt(2)) & all(rt_vals>xlimhist_rt(1)))
            
            histogram(rt_vals,[xlimhist_rt(1):0.05:xlimhist_rt(2)],'FaceColor',task_colors(tt,:),'EdgeColor','none','FaceAlpha',0.4);

            xlim(xlimhist_rt);
            ylim(ylimhist_rt);
            line([0,0],get(gca,'YLim'),'Color','k')
            line([meanrt,meanrt],get(gca,'YLim'),'Color',task_colors(xx,:),'LineWidth',3)
            line([meanrt-stdev,meanrt-stdev],get(gca,'YLim'),'Color',task_colors(xx,:))
            line([meanrt+stdev,meanrt+stdev],get(gca,'YLim'),'Color',task_colors(xx,:))
            set(gcf,'Color','w');
            xlabel('RT (s)');
            ylabel('Number of trials');
            set(gca,'FontSize',12);
      
    end
    suptitle(sprintf('Response time, all MRI subjects (n=%d)',nSubj))

    set(gcf,'Color','w');
    set(gca,'FontSize',12);
    size = get(0,'ScreenSize');
    set(gcf,'Position',[200,200,800,1400])
    set(gcf,'Color','w')
    
end
    
    
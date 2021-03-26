%% analyze behavior of all subjects during 3 tasks

clear
close all

currdir = pwd;
ndirsup=1;
filesepinds = find(currdir==filesep);
root = currdir(1:filesepinds(end-ndirsup+1));
addpath('/usr/local/serenceslab/maggie/mFiles/Classifiers/');

%% LOADING DATA FROM VERSION 2

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

acc_by_task_V2 = nan(nSubj, nTasks, 8);
acc_by_time_V2 = nan(nSubj, 2, 6, 2);
legendlabs = [];

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
nRespRight_V2 = zeros(nSubj,2,nDistances);
nRespLeft_V2 = zeros(nSubj,2,nDistances);
nRespFail_V2 = zeros(nSubj,2,nDistances);

% nRespGrid = zeros(nSubj, nTasks, numel(all_pts),numel(all_pts), 5);

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
                
                acc_by_task_V2(si,curr_task,task_counts(curr_task)) = acc;
                acc_by_time_V2(si,se,xx,rr) = acc;
                
                
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
                            nRespRight_V2(si,curr_task,ind) = nRespRight_V2(si,curr_task,ind) + 1;
                        elseif resp_unmapped==1
                            nRespLeft_V2(si,curr_task,ind) = nRespLeft_V2(si,curr_task,ind) + 1;
                        else
                            nRespFail_V2(si,curr_task,ind) = nRespFail_V2(si,curr_task,ind) + 1;
                        end
                        
%                         nansum(nRespFail(si,curr_task,:))+nansum(nRespRight(si,curr_task,:))+nansum(nRespLeft(si,curr_task,:))
                        
                    end
                    assert(nansum(nRespFail_V2(si,curr_task,:))+nansum(nRespRight_V2(si,curr_task,:))+nansum(nRespLeft_V2(si,curr_task,:))==nTrials*task_counts(curr_task))

                end
               
            end
        end
    end
    
    %% check a few things

    for tt=1:nTasks
        if tt<3
            assert(nansum(nRespFail_V2(si,tt,:))+nansum(nRespRight_V2(si,tt,:))+nansum(nRespLeft_V2(si,tt,:))==nTrials*task_counts(tt))
        end
      
    end
end

%% LOADING DATA FROM VERSION 4

% first list ALL subjects 
sublist = [1:8];

% now list any that I want to remove
sub2remove = [];
inds2remove = ismember(sublist,sub2remove);
sublist(inds2remove) = [];

% how many are left?
nSubj = length(sublist);

% list the orders in which things will happen across entire experiment
% [nParts x nSess]
task_list_odd = [1,2;2,1;3,3]; % which category boundary?
map_list_odd = [1,2;1,2;1,2];   % which response mapping?
position_list_odd = [1,2;1,2;1,2];  % which positions will we show prototypes on the screen?
% even num subject, start with second sess (otherwise start with
% the first session)
task_list_even = task_list_odd(:,[2,1]);
map_list_even = map_list_odd(:,[2,1]);
position_list_even = position_list_odd(:,[2,1]);

% names of the tasks
tasklist = {'Linear (1)','Linear (2)','Checker'};
nTasks = length(tasklist);

nRunsPerPart = 3;
nPartsPerTask = 1;  % 2 mappings per task
nPartsTotal  = nTasks*nPartsPerTask;
nSess = 2;
% acc_by_task_V4 = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_by_task_V4 = nan(nSubj, nTasks, 8);
acc_by_run = nan(nSubj, nSess, nPartsTotal, nRunsPerPart);

legendlabs = [];

plotRTs = 1;
plotAccByTime = 1;
plot2DChoice=1;
plot1DChoice=1;

%% quality control settings...
% what percent of trials can be timeout before i toss the whole block?
% if they fail to respond to 25 percent, probably they were confused and/or
% distracted and/or unmotivated.
timeout_pct_cutoff = 0.25;

%% create point grid in feature space

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  

distances_from_bound = [-0.8:0.1:-0.1, 0.1:0.1:0.8];
nDistances = numel(distances_from_bound);

% create psychometric curves - how likely is it that they report to the
% right/left of center boundary? [subj,task,position along relevant axis]
% for linear tasks only
nRespRight_V4 = zeros(nSubj,2,nDistances);
nRespLeft_V4 = zeros(nSubj,2,nDistances);
nRespFail_V4 = zeros(nSubj,2,nDistances);

nRespGrid = zeros(nSubj, nTasks, numel(all_pts),numel(all_pts), 3);

possible_responses =[1,2];
%% load data

actual_cb = zeros(nSubj,1);
for si = 1:nSubj

    subinit = sprintf('S%02.f', sublist(si));
    
    % look at all data in the folder, figure out the dates of the two
    % sessions
    alldat = dir(fullfile(root, 'Data_V4',sprintf('%s_MainTask_*.mat',subinit)));
    names = {alldat.name};
    dates = cellfun(@(x)x(end-9:end-4),names,'UniformOutput',false);
    sess_dates = unique(dates);
    [~,order] = sort(str2double(sess_dates),'ascend');
    sess_dates = sess_dates(order);
    nSess = numel(sess_dates);
    if nSess~=2
        error('%s: has %d sessions\n',subinit,nSess)
    end
    SubNum=sprintf('%02d',sublist(si));
    % figure out which counter-balance condition this subject was in
    if ~mod(SubNum,2)
        actual_cb(si) = 2;
        %     if mod(si,2)
        fprintf('%d: even\n',si);
        task_list = task_list_even;
        map_list = map_list_even;
        position_list = position_list_even;
    else
        actual_cb(si) = 1;
        fprintf('%d: odd\n',si);
        task_list = task_list_odd;
        map_list = map_list_odd;
        position_list = position_list_odd;        
    end
    
    legendlabs{si} = sprintf('%s',subinit);
    
    % counting tasks - combining both parts of each task together 
    task_counts = zeros(nTasks,1);
    
    allrts = [];
    
    % loop over sessions
    for se = 1:nSess
        
        % looping over the 6 parts of this session
        for xx = 1:nPartsTotal
            
            filename = fullfile(root, 'Data_V4',sprintf('%s_MainTask_sess%d_part%d_%s.mat',subinit,se,xx,sess_dates{se}));
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
                
                acc_by_task_V4(si,curr_task,task_counts(curr_task)) = acc;
                acc_by_run(si,se,xx,rr) = acc;
                
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
                        dist = round(p.difficulty(tr),1);
                        
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
                            nRespRight_V4(si,curr_task,ind) = nRespRight_V4(si,curr_task,ind) + 1;
                        elseif resp_unmapped==2
                            nRespLeft_V4(si,curr_task,ind) = nRespLeft_V4(si,curr_task,ind) + 1;
                        else
                            nRespFail_V4(si,curr_task,ind) = nRespFail_V4(si,curr_task,ind) + 1;
                        end
                        
%                         nansum(nRespFail(si,curr_task,:))+nansum(nRespRight(si,curr_task,:))+nansum(nRespLeft(si,curr_task,:))
                        
                    end
                    assert(nansum(nRespFail_V4(si,curr_task,:))+nansum(nRespRight_V4(si,curr_task,:))+nansum(nRespLeft_V4(si,curr_task,:))==nTrials*task_counts(curr_task))

                end
                
                
            end
            
            
        end
        
    end
   
    %% check a few things

    for tt=1:nTasks
        if tt<3
            assert(nansum(nRespFail_V4(si,tt,:))+nansum(nRespRight_V4(si,tt,:))+nansum(nRespLeft_V4(si,tt,:))==nTrials*task_counts(tt))
        end
       
    end
end

%% COMBINING DATA FROM BOTH VERSIONS OF BINARY TASKS

acc_by_task = cat(1, acc_by_task_V2, acc_by_task_V4);

nRespFail = cat(1, nRespFail_V2, nRespFail_V4);
nRespLeft = cat(1, nRespLeft_V2, nRespLeft_V4);
nRespRight = cat(1, nRespRight_V2, nRespRight_V4);

nSubj1 = size(acc_by_task_V2,1);
nSubj2 = size(acc_by_task_V4,1);

nSubj = size(acc_by_task,1);

%% plot accuracy at each task
figure();hold all;
sub_colors = viridis(nSubj);
vals = nanmean(acc_by_task(:,1:2,:), 3);

legendlabs=[];
for si=1:nSubj
    if si<nSubj1
        legendlabs{si} = sprintf('S%02d, V2',si);
    else
        legendlabs{si} = sprintf('S%02d, V4',si-nSubj1+1);
    end
    plot(1:2, vals(si,1:2),'-o','Color',sub_colors(si,:));

end

meanvals = squeeze(nanmean(vals,1));
sevals = squeeze(nanstd(vals,[],1))./sqrt(nSubj);
errorbar(1:2, meanvals,sevals,'Color','k');
plot([0,4],[0.5, 0.5],'Color','k','LineStyle','--');
xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Accuracy, averaged over all runs'))

set(gca,'XLim',[0,3], 'YLim',[0.0, 1]);
set(gca,'XTick',1:2,'XTickLabel',tasklist(1:2),'XTickLabelRotation',45)

set(gcf,'Color','w');

legend(legendlabs,'Location','EastOutside')



%% plot psychometric curves (only for the linear tasks)

cols = viridis(nSubj);

% legendlabs=[];
if plot1DChoice
    
    for tt=1:2
        
        figure();hold all;
        props = squeeze(nRespRight(:,tt,:)./(nRespLeft(:,tt,:)+nRespRight(:,tt,:)));
        
        for si=1:nSubj
            
%             if si<nSubj1
%                 legendlabs{si} = sprintf('S%02d, V2',si);
%             else
%                 legendlabs{si} = sprintf('S%02d, V4',si-nSubj1+1);
%             end
            plot(distances_from_bound,squeeze(props(si,:)),'-','Color',[0.8, 0.8, 0.8])
%             plot(distances_from_bound,squeeze(props(si,:)),'-','Color',cols(si,:))
            
        end
        
        if nSubj>1
            meanvals = squeeze(nanmean(props,1));
            sevals = squeeze(nanstd(props,[],1))./sqrt(nSubj);
        else
            meanvals = props;
            sevals = [];
        end
        errorbar(distances_from_bound, meanvals,sevals,'Color','k');
        
        
        xlabel(sprintf('distance from center'))
        ylabel(sprintf('Proportion responded right'));
        title(sprintf('Psychometric curve, mean over runs\n %s task',tasklist{tt}))
        xlim([-1,1])
        ylim([0,1])
        plot([0,0], get(gca,'YLim'),'k')
        %
        plot(get(gca,'XLim'), [0.5,0.5],'k')
        set(gcf,'Color','w')
    end
    
    
end
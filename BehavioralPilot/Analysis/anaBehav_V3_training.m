%% analyze behavior of all subjects during 3 tasks

clear
close all

root = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/BehavioralPilot/';
addpath('/usr/local/serenceslab/maggie/mFiles/Classifiers/');

% first list ALL subjects 
sublist = [1:13];
% sublist = [];

% now list any that I want to remove
sub2remove = [];
inds2remove = ismember(sublist,sub2remove);
sublist(inds2remove) = [];

% how many are left?
nSubj = length(sublist);

% list the orders in which things will happen across entire experiment
% [nParts x nSess]
task_list_odd = [1,2;1,2;2,1;2,1;3,3;3,3]; % which category boundary?
map_list_odd = [1,2;2,1;1,2;2,1;1,2;2,1];   % which response mapping?
position_list_odd = [1,2;1,2;1,2;1,2;1,2;1,2];  % which positions will we show prototypes on the screen?
% even num subject, start with second sess (otherwise start with
% the first session)
task_list_even = task_list_odd(:,[2,1]);
map_list_even = map_list_odd(:,[2,1]);
position_list_even = position_list_odd(:,[2,1]);

% names of the tasks
tasklist = {'Linear (1)','Linear (2)','Checker'};
nTasks = length(tasklist);


nRunsPerPart = 2;
nPartsPerTask = 2;  % 2 mappings per task
nPartsTotal  = nTasks*nPartsPerTask;
nSess = 2;
acc_by_task = nan(nSubj, nTasks, nRunsPerPart*nPartsPerTask*nSess);
acc_by_run = nan(nSubj, nSess, nPartsTotal, nRunsPerPart);

legendlabs = [];

plotRTs = 0;
plotAccByTime = 1;
plot2DChoice=1;
plot1DChoice=1;
%% create point grid in feature space

% first I'm defining all possible images that we can use in this task.
start = 0.2;    % min value along each axis
stop = 4.8; % max value along each axis  
step = 0.1;
center = (stop-start)./2+start;

all_pts = start:step:stop;  
% [gridx,gridy] = meshgrid(all_pts,all_pts);
% all_grid_points = [gridx(:),gridy(:)];

    % now taking out images at exactly the prototype locations so that we
%     % never use these during task
%     nsteps_main = 4;
%     main_pts = round(linspace(start,stop, nsteps_main),1);
%     proto_pts = [round(mean(main_pts(1:2)),1), round(mean(main_pts(3:4)),1)];
%     proto_coords = [proto_pts(2), proto_pts(2); proto_pts(1), proto_pts(2); proto_pts(1), proto_pts(1); proto_pts(2), proto_pts(1)];
%     proto_inds = find(ismember(all_grid_points, proto_coords, 'rows'));
%     assert(numel(proto_inds)==4)
%     all_grid_points(proto_inds,:) = [];
%   
%     % Also taking out any images along quadrant boundaries because these
%     % are ambiguous 
%     bound_inds = find(any(all_grid_points==center,2));
%     all_grid_points(bound_inds,:) = [];
%   
%     % now define which quadrant each image lies in
%     all_quadrant = zeros(size(all_grid_points,1),1);
%     all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)>center) = 1;
%     all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)>center) = 2;
%     all_quadrant(all_grid_points(:,1)<center & all_grid_points(:,2)<center) = 3;
%     all_quadrant(all_grid_points(:,1)>center & all_grid_points(:,2)<center) = 4;

distances_from_bound = [-0.8:0.1:-0.1, 0.1:0.1:0.8];
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
    
    % look at all data in the folder, figure out the dates of the two
    % sessions
    alldat = dir(fullfile(root, 'Data_V3',sprintf('%s_ShapeTask_*TRAINING*.mat',subinit)));
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
            
            filename = fullfile(root, 'Data_V3',sprintf('%s_ShapeTask_sess%d_part%d_TRAINING_%s.mat',subinit,se,xx,sess_dates{se}));
            load(filename);
            p = TheData(1).p;
            
            if se==1 && xx==1
                if ~mod(p.SubNum,2)
                    fprintf('%s: even\n',subinit);
                else
                    fprintf('%s: odd\n',subinit);
                end
            end
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
                    fprintf('Subject %s is missing run %d for sess %d part %d\n',p.SubNum,rr,se,xx);
                    continue
                end
                
                p = TheData(rr).p;
                data = TheData(rr).data;
                t = TheData(rr).t;
                
                % keep track of how many runs of the current task i have 
                task_counts(curr_task) = task_counts(curr_task) + 1;
%                 assert(task_counts(curr_task)==rr)
                
                timeout = data.Response==0 | ~ismember(data.Response,possible_responses);
                
                if sum(timeout)>5
                    % print a warning if there were lots of timeouts
                    fprintf('%s sess %d task %d map %d run %d: %d timeout trials\n',subinit,se,curr_task,curr_map,rr,sum(timeout));
                end
                
                rts = t.RespTimeFromOnset;
                rts = rts(~timeout);
                allrts = [allrts; rts];
                
%                 resp = p.category(~timeout)';
                resp = data.Response(~timeout);
                correct_resp = p.category(~timeout)';
                acc = mean(resp==correct_resp);
                
                acc_by_task(si,curr_task,task_counts(curr_task)) = acc;
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

vals = nanmean(acc_by_task, 3);

for si=1:nSubj
   
    plot(1:3, vals(si,:),'-o');

end

meanvals = squeeze(nanmean(vals,1));
sevals = squeeze(nanstd(vals,[],1))./sqrt(nSubj);
errorbar(1:3, meanvals,sevals,'Color','k');

xlabel(sprintf('Task'))
ylabel(sprintf('Accuracy'));
title(sprintf('Training Accuracy, averaged over runs'))

set(gca,'XLim',[0,4], 'YLim',[0.0, 1]);
set(gca,'XTick',1:3,'XTickLabel',tasklist,'XTickLabelRotation',45)

set(gcf,'Color','w');

legend(legendlabs,'Location','EastOutside')


%% plot accuracy in order of when blocks were done
if plotAccByTime
for se = 1:2
    
    for cb = 1:2

        figure();hold all;

        sub2use = find(actual_cb'==cb);
        
        for si=sub2use

            dat = squeeze(acc_by_run(si,se,:,:))';
            dat = dat(:);
            plot(1:numel(dat), dat,'-o');

        end

        xlabel(sprintf('Task'))
        ylabel(sprintf('Accuracy'));
        title(sprintf('Session %d, cb %d\nTraining task accuracy versus run number',se,cb))

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
if plot2DChoice
    
allax=  [];
tick_pts = [1,24,47];

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
%             plot(1, 2.5, 'ro')
        colorbar()
        title(sprintf('probability of resp %d',rr))
        axis square 
        set(gca,'XTick',tick_pts,'XTickLabel',all_pts(tick_pts));
        set(gca,'YTick',tick_pts,'YTickLabel',all_pts(tick_pts),'YDir','normal');
        xlabel('feature 1');
        ylabel('feature 2');
    end
    
    suptitle(sprintf('%s TRAINING task',tasklist{tt}))
    
    set(gcf,'Color','w')
end

match_clim(allax)

end
%% plot psychometric curves (only for the linear tasks)

cols = viridis(nSubj);

if plot1DChoice
    
    for tt=1:2
        
        figure();hold all;
        props = squeeze(nRespRight(:,tt,:)./(nRespLeft(:,tt,:)+nRespRight(:,tt,:)));
        
        for si=1:nSubj
            
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
        title(sprintf('Psychometric curve, mean over runs\n %s TRAINING task',tasklist{tt}))
        xlim([-1,1])
        ylim([0,1])
        plot([0,0], get(gca,'YLim'),'k')
        %
        plot(get(gca,'XLim'), [0.5,0.5],'k')
        set(gcf,'Color','w')
    end
    
    
end
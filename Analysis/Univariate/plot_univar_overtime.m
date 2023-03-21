%% Plot average signal in each ROI over time. 
% similar to deconvolution but instead of deconvolving, just stack epoched
% trials and take average. 

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1));

%% set up some parameters

zscore_each_run = 1;

nTasks = 3;
task_names = {'Linear (1)','Linear (2)','Checker','Repeat Detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);


nMainRunsTotal = nTasks*2*2*3;
nRepRunsTotal = 4*3;
nTrialsPerRun = 48;

ROI_names = {'V1','V2','V3','V3AB','hV4','IPS0','IPS1','IPS2','IPS3','LO1','LO2',...
    'IFS', 'AI-FO', 'iPCS', 'sPCS','sIPS','ACC-preSMA'};
nROIs = length(ROI_names);

plot_order_all = 1:length(ROI_names);
plot_names = ROI_names(plot_order_all);

trDur = .8; % actual duration of the TR, for plotting, in seconds...
nTRs_main = 327-16;
nTRs_rep = 329-16;

% events to plot as vertical lines
evt_times_plot = [0,1];
nTRs_out = 14;
% define timing axis
tax = trDur*(0:nTRs_out-1); 

%%

allsub_bold_mean = nan(nSubj, length(plot_order_all), nTasks+1, nTRs_out);
allsub_bold_sem = nan(nSubj, length(plot_order_all), nTasks+1, nTRs_out);

for ss=1:nSubj
   
    substr = sprintf('S%02d',ss);
   
    fn2load = fullfile(root,'Samples',sprintf('SampleFile_%s.mat',substr));
    load(fn2load, 'samplesMain','samplesRep','ROIs','all_vox_concat');
     
    %% load the timing file (made in GetEventTiming.m)
    
    fn = fullfile(root,'Samples',sprintf('TimingFile_%s.mat',substr));
    if ~exist(fn, 'file')
        error('need to make timing file first, run GetEventTiming.m')
    end
    fprintf('Loading event timing file\n')
    load(fn)
    
    for vv = 1:nROIs % for all visual areas I want to look at
        %% pull out the data from each ROI
        % want both hemispheres
        [rowind1,colind1] = find(strcmp(reshape({ROIs.name},2,nROIs),sprintf('lh_%s',ROI_names{vv})));
        [rowind2,colind2] = find(strcmp(reshape({ROIs.name},2,nROIs),sprintf('rh_%s',ROI_names{vv})));
        col_inds = [colind1,colind2]; % column is the region
        row_inds = [rowind1,rowind2];   % row is the hemisphere
        
        if vv<12
            % make sure each early visual area is defined for both
            % hemispheres.
            assert(numel(col_inds)==2)
        end
        
        mainDat=[];
        repDat = [];
        for ii=1:length(col_inds)
            name = ROIs(row_inds(ii),col_inds(ii)).name;
            if ~isempty(ROIs(row_inds(ii),col_inds(ii)).voxel_inds)
                % jj gives indices into the all_vox_concat array
                [~,jj]=intersect(all_vox_concat, ROIs(row_inds(ii),col_inds(ii)).voxel_inds);
                mainDat = [mainDat, samplesMain(:,jj)];
                repDat = [repDat, samplesRep(:,jj)];
            end
        end
        nVox = size(mainDat,2);
        assert(size(repDat,2)==nVox);

        if nVox==0
            fprintf('no voxels in area %s!\n',ROI_names{vv});
            continue
        end

        fprintf('processing area %s, %d voxels\n', ROI_names{vv}, nVox);
        
        %%
        if zscore_each_run
            nRunsMain = size(mainDat,1)/nTRs_main;
            assert(nRunsMain==nMainRunsTotal);
            if mod(nRunsMain,1)~=0
                error('something bad happened here with mainDat run length')
            end
            for ii=1:nRunsMain
                mainDat(ii*nTRs_main-nTRs_main+1:ii*nTRs_main,:) = zscore(mainDat(ii*nTRs_main-nTRs_main+1:ii*nTRs_main, :),1);
            end

            assert(numel(unique(main.RunLabels))==nRunsMain); 
            
            nRunsRep = size(repDat,1)/nTRs_rep;
            assert(nRunsRep==nRepRunsTotal);
            if mod(nRunsRep,1)~=0
                error('something bad happened here with repDat run length')
            end
            for ii=1:nRunsRep
                repDat(ii*nTRs_rep-nTRs_rep+1:ii*nTRs_rep,:) = zscore(repDat(ii*nTRs_rep-nTRs_rep+1:ii*nTRs_rep, :),1);
            end

            assert(numel(unique(rep.RunLabels))==nRunsRep); 
        end
        
        %% label the data
        % 1 = stim on, 0 = stim off
        event_diff_reshaped = reshape(diff([0;main.EventLabels]),nTRs_main,length(main.EventLabels)/nTRs_main);

        % now find the actual onset of each stimulus - switch from 0.2 to 1
        % (or 0 to 1)
        trial_onset_bool = event_diff_reshaped==1;
        trial_onset_bool = trial_onset_bool(:);
        trial_onset_num = find(trial_onset_bool);
                
        nTrialsMain = nMainRunsTotal*nTrialsPerRun;
        assert(numel(trial_onset_num)==nTrialsMain);
        
        task_labels = main.BoundLabels(trial_onset_num);
        %% gather data each trial
        
        % TR-by-TR data        
        maindat_by_TR = nan(nTrialsMain, nTRs_out, nVox);
        
        triCnt = 0; % counter across "good" trials where the entire desired avg window is available.
    
        for rr=1:nMainRunsTotal

            runInd = rr*nTRs_main-nTRs_main+1:rr*nTRs_main;
            assert(all(find(main.RunLabels==rr)==runInd'))
            curDat = mainDat(runInd,:);      % data from current run

            % get numeric indices for each event in this run
            these_targ_onsets = find(trial_onset_bool(runInd));                   

            assert(numel(these_targ_onsets)==nTrialsPerRun);
            
            for tt=1:numel(these_targ_onsets)
                triCnt = triCnt + 1;  % increment counter over good trials
                % also collecting data at each timept
                for tr = 1:nTRs_out
                    if these_targ_onsets(tt)+tr-1<=nTRs_main 
                        maindat_by_TR(triCnt,tr,:) = curDat(these_targ_onsets(tt)+tr-1,:);
                    else
                        error('TR %d is past the end of your run!!', tr)
                    end
                end
                
            end
        end

        assert(triCnt==nTrialsMain)
        assert(~sum(isnan(maindat_by_TR(:))))
        
        %% same thing for repeat detection task
        % 1 = stim on, 0 = stim off
        event_diff_reshaped = reshape(diff([0;rep.EventLabels]),nTRs_rep,length(rep.EventLabels)/nTRs_rep);

        % now find the actual onset of each stimulus - switch from 0.2 to 1
        % (or 0 to 1)
        trial_onset_bool = event_diff_reshaped==1;
        trial_onset_bool = trial_onset_bool(:);
        trial_onset_num = find(trial_onset_bool);
                
        nTrialsRep = nRepRunsTotal*nTrialsPerRun;
        assert(numel(trial_onset_num)==nTrialsRep);
        
        %% gather data each trial
        
        % TR-by-TR data        
        repdat_by_TR = nan(nTrialsRep, nTRs_out, nVox);
        
        triCnt = 0; % counter across "good" trials where the entire desired avg window is available.
    
        for rr=1:nRepRunsTotal

            runInd = rr*nTRs_rep-nTRs_rep+1:rr*nTRs_rep;
            assert(all(find(rep.RunLabels==rr)==runInd'))
            curDat = repDat(runInd,:);      % data from current run

            % get numeric indices for each event in this run
            these_targ_onsets = find(trial_onset_bool(runInd));                   
          
            
            assert(numel(these_targ_onsets)==nTrialsPerRun);
            
            for tt=1:numel(these_targ_onsets)
                triCnt = triCnt + 1;  % increment counter over good trials
                % also collecting data at each timept
                for tr = 1:nTRs_out
                    if these_targ_onsets(tt)+tr-1<=nTRs_rep 
                        repdat_by_TR(triCnt,tr,:) = curDat(these_targ_onsets(tt)+tr-1,:);
                    else
                        error('TR %d is past the end of your run!!', tr)
                    end
                end
                
            end
        end

        assert(triCnt==nTrialsRep)
        assert(~sum(isnan(repdat_by_TR(:))))
       
        %% now group into tasks and average
        
        for tt = 1:nTasks
            mean_over_trials = squeeze(mean(maindat_by_TR(task_labels==tt,:,:), 1));
            
            allsub_bold_mean(ss,vv,tt,:) = mean(mean_over_trials, 2);
            allsub_bold_sem(ss,vv,tt,:) = std(mean_over_trials, [], 2)/sqrt(nVox);
        end
        
        tt=nTasks+1;
        mean_over_trials = squeeze(mean(repdat_by_TR, 1));            
        allsub_bold_mean(ss,vv,tt,:) = mean(mean_over_trials, 2);
        allsub_bold_sem(ss,vv,tt,:) = std(mean_over_trials, [], 2)/sqrt(nVox);
        
    end
end


%%
ylims = [-0.25, 0.25];
lw=1;
nplots = ceil(sqrt(numel(plot_order_all)));
figure;hold all;

for vi = 1:numel(plot_order_all)
    subplot(nplots,ceil(numel(plot_order_all)/nplots),vi);hold all;
    lh=[];lx=0;ll=[];
    for tt=1:nTasks+1

        % take contra - ipsi difference
        vals = allsub_bold_mean(:,plot_order_all(vi),tt,:);

        if nSubj==1
            meanvals =squeeze(vals);
            semvals = zeros(size(meanvals));
%             semvals = squeeze(allsub_bold_sem(:,plot_order_all(vi),tt,:));
        else
            meanvals = squeeze(nanmean(vals,1));
            semvals = squeeze(nanstd(vals,[],1)./sqrt(sum(~isnan(vals),1)));
        end

        lh=[lh,plot(tax,meanvals,'-','Color',task_colors(tt,:),'LineWidth',lw)];
        bandedError_MMH(tax, meanvals',semvals', task_colors(tt,:), 0.5);
%                 lh=[lh, errorbar(tax, meanvals, semvals,'Color',col_conds(cc,:,bb),'LineWidth',lw)];
        lx=lx+1;
        ll{lx} = sprintf('%s',task_names{tt});

    end


    set(gca, 'FontSize', 12, 'XLim',[0 max(tax)],'XTick',[0,max(tax)/2, max(tax)])

    if zscore_each_run
            ylim(ylims)
            h=plot(get(gca,'XLim'),[0,0],'--','Color',[0.7, 0.7, 0.7]);
            uistack(h,'bottom')
            for ee = 1:length(evt_times_plot)
                h=plot([evt_times_plot(ee),evt_times_plot(ee)],ylims,'-','Color',[0.8, 0.8, 0.8]);
                uistack(h,'bottom')
            end
            h=fill([evt_times_plot(1),evt_times_plot(2), evt_times_plot(2),evt_times_plot(1)], repelem(ylims,2),[0.8, 0.8, 0.8],'EdgeColor','none');
            uistack(h,'bottom')
    
    end
%     if vi==1
%         
%     end

    if vi==numel(plot_order_all)
        legend(lh, ll);
        xlabel('Time(s)');
        ylabel('BOLD resp');
    end


    title(sprintf('%s',plot_names{vi}));

end

if zscore_each_run
    suptitle('Mean Z-scored BOLD, all trials each condition')
else
    suptitle('Mean BOLD without z-scoring, all trials each condition')
end
set(gcf,'Color','w')
set(gcf,'Position',[200,200,1800,1200]);


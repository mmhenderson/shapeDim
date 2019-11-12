%% analyze behavior on pilot task

clear
close all

% set inputs
subnum = 1;

root = pwd;
filesepinds = find(root==filesep);
root = root(1:filesepinds(end-1));

nRunsPerSet = 8;
nSets = 2;
acc = zeros(nRunsPerSet, nSets);
dprime = zeros(nRunsPerSet, nSets);

nIms = 16;
confmat = zeros(nIms, nIms, 3, nSets);

% set paths (end all with filesep!)
experiment_path = fullfile(root,'Pilot1');
% get data file from each session (they are not concatenated here)
loc_file = dir(fullfile(experiment_path,'DataBehavior',sprintf('S%02d',subnum),'Session*','*OneBackPilotTask*'));
load(fullfile(loc_file(1).folder,loc_file(1).name));
nTrialsEach = TheData(1).p.nTrials;

run_set_count = zeros(nSets,1);

for run = 1:length(TheData)

    ss = TheData(run).p.stimset-2;
    run_set_count(ss) = run_set_count(ss)+1;
    
    for n = 2:TheData(run).p.nTrials
   
        prev_im = TheData(run).p.imlist(n-1);
        curr_im = TheData(run).p.imlist(n);
        resp_ind = TheData(run).data.Response(n);
        if isnan(resp_ind) || resp_ind==0
            resp_ind=3;
        end
        confmat(prev_im, curr_im, resp_ind, ss) = confmat(prev_im, curr_im, resp_ind, ss) + 1;
       
    end
    
    resp = TheData(run).data.Response;
    correct_resp = TheData(run).p.correct_resp;
    inds2use = find(TheData(run).data.Response~=0 & ~isnan(TheData(run).data.Response));
    
    thisacc = mean(resp(inds2use)==correct_resp(inds2use));
    acc(run_set_count(ss), ss) = thisacc;
    
    if any(correct_resp(inds2use)==1)
        thisd = get_dprime(resp(inds2use),correct_resp(inds2use),unique(correct_resp));
        dprime(run_set_count(ss),ss) = thisd;
    else
        dprime(run_set_count(ss),ss) = nan;
    end
end %run loop

%%
colormap = plasma(3);
figure;hold all;
legendlabs = [];
for ss=1:size(acc,2)
    plot(1:nRunsPerSet, acc(:,ss),'-o','MarkerEdgeColor',colormap(ss,:),'MarkerFaceColor',colormap(ss,:));
    legendlabs = [legendlabs, {sprintf('set %d',ss+2)}];
end
% set(gca,'XTick',[1,2],'XTickLabel',{'Set 3','Set 4'});
xlabel('Run number');
xlim([0,9]);
ylabel('Accuracy');
ylim([0.6, 1.1]);
legend(legendlabs, 'Location','EastOutside');
title('accuracy on one-back repeat detection task over time');
%%
colormap = plasma(3);
figure;hold all;
legendlabs = [];
for ss=1:size(acc,2)
    plot(1:nRunsPerSet, dprime(:,ss),'-o','MarkerEdgeColor',colormap(ss,:),'MarkerFaceColor',colormap(ss,:));
    legendlabs = [legendlabs, {sprintf('set %d',ss+2)}];
end
% set(gca,'XTick',[1,2],'XTickLabel',{'Set 3','Set 4'});
xlabel('Run number');
xlim([0,9]);
ylabel('d-prime');
% ylim([0.6, 1.1]);
legend(legendlabs, 'Location','EastOutside');
title('d-prime on one-back repeat detection task over time');

%%

prop_match = squeeze(confmat(:,:,1,:))./squeeze(sum(confmat, 3));
not_tested = squeeze(sum(confmat,3))==0;
prop_match(not_tested) = nan;

tick_pts = 1:nIms;

for ss=1:nSets
   
    figure;hold all;
    sanePColor(1:nIms,1:nIms,prop_match(:,:,ss));
%             plot(1, 2.5, 'ro')
    colorbar()
    title(sprintf('probability of resp match, set %d',ss+2))
    axis square 
    set(gca,'XTick',tick_pts);
    set(gca,'YTick',tick_pts,'YDir','normal');
    xlabel('current image');
    ylabel('prev image');
end
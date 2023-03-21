% measure overall decodability of shapes

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1));
% savepath = fullfile(root, 'Analysis','Decoding_results','BinaryAcc_SVM_all.mat');
% savepath = fullfile(root, 'Analysis','Decoding_results','Decode_task_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Task: Linear (1)','Task: Linear (2)','Task: Checker','Repeat detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);

nRunsEach = 2*2*3;

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
%         all_acc = zeros(nSubj, nVOIs, nPairs);
        avg_signal = zeros(nSubj, nVOIs,nTasks+1, nRunsEach);

    end

    for vv=1:nVOIs
        
        dat = mainSig(vv).dat_avg;
        
        for tt=1:nTasks

            unruns = unique(mainSig(vv).RunLabels(mainSig(vv).BoundLabels==tt));
            assert(numel(unruns)==nRunsEach);
            
            for rr=1:nRunsEach
                
                inds2use = mainSig(vv).BoundLabels==tt & mainSig(vv).RunLabels==unruns(rr);
                assert(sum(inds2use)==48);
                
                avg_signal(ss,vv,tt,rr) = mean(mean(dat(inds2use,:),1),2);
                
            end
        end
        
        dat = repSig(vv).dat_avg;
        tt=nTasks+1;
        unruns = unique(repSig(vv).RunLabels);
        assert(numel(unruns)==nRunsEach);

        for rr=1:nRunsEach

            inds2use = repSig(vv).RunLabels==unruns(rr);
            assert(sum(inds2use)==48);

            avg_signal(ss,vv,tt,rr) = mean(mean(dat(inds2use,:),1),2);
            
        end
    end
end


%% plot decoding along each axis separately
v2plot = [1:5,10,11,6:9,12:17];


mean_acc = squeeze(mean(mean(avg_signal(:,v2plot,:,:),4),1));
se_acc = squeeze(mean(std(avg_signal(:,v2plot,:,:),[],4),1))./sqrt(nRunsEach);

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
% se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],0,[],ROI_names(v2plot),task_names,'Avg z-scored BOLD',[],task_colors)
set(gcf,'Position',[200,200,1000,400])

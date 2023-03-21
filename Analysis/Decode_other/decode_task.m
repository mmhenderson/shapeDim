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
savepath = fullfile(root, 'Analysis','Decoding_results','Decode_task_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Task: Linear (1)','Task: Linear (2)','Task: Checker'};
task_pairs = combnk(1:nTasks, 2);
nPairs = size(task_pairs,1);

nRunsTotal=nTasks*2*2*3;
nTrialsPerRun = 48;

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nPairs);
    end

    for vv=1:nVOIs

        for pp = 1:nPairs
            
            tasklabs = mainSig(vv).BoundLabels;
%             maplabs = mainSig(vv).MapLabels;
%             inds2use = maplabs==1 & (tasklabs==task_pairs(pp,1) | tasklabs==task_pairs(pp,2));
            inds2use = tasklabs==task_pairs(pp,1) | tasklabs==task_pairs(pp,2);

            dat = mainSig(vv).dat_avg(inds2use,:);
    
            tasklabs = mainSig(vv).BoundLabels(inds2use);
            
            un = unique(tasklabs);
            num_each = sum(repmat(tasklabs,1,numel(un))==repmat(un',numel(tasklabs),1),1);
            assert(all(num_each==num_each(1)))
         
            cvlabs = mainSig(vv).SessLabels(inds2use);
         
            nCV = numel(unique(cvlabs));
            
            predlabs = nan(size(tasklabs));
            
            for cv=1:nCV

                trninds = cvlabs~=cv;
                tstinds = cvlabs==cv;
                
                trnlabs = tasklabs(trninds);

                un = unique(trnlabs);
                num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
                assert(all(num_each==num_each(1)))
                
%                 label = classify(dat_thispair(tstinds,:),dat_thispair(trninds,:),pairlabs,'diaglinear');
                
                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
    %             m = fitcsvm(dat(trninds,:), tasklabs(trninds));
    %             label = predict(m, dat(tstinds,:));
                predlabs(tstinds) = label;

            end
            assert(~any(isnan(predlabs)))

            acc = mean(predlabs==tasklabs);
            fprintf('S%s, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},acc*100);
            all_acc(ss,vv,pp) = acc;

        end

        

    end
    
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)
%% plot decoding along each axis separately
v2plot = [1:5,10,11,6:9,12:17];


mean_acc = squeeze(mean(mean(all_acc(:,v2plot,:),3),1))';

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],1/2,[0.1, 1],ROI_names(v2plot),[],'Accuracy',[],[])
set(gcf,'Position',[200,200,1000,400])

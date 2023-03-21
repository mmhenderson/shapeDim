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
savepath = fullfile(root, 'Analysis','Decoding_results','Decode_respfinger_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Task: Linear (1)','Task: Linear (2)','Task: Checker','Repeat Detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);


nRunsTotal=nTasks*2*2*3;
nTrialsPerRun = 48;

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks+1);
    end

    for vv=1:nVOIs
        
        for tt=1:nTasks
 
            inds2use = mainSig(vv).BoundLabels==tt;
%   
            dat = mainSig(vv).dat_avg(inds2use,:);
            dat = dat - repmat(mean(dat,2),1,size(dat,2));
  
            fingerlabs = mainSig(vv).CatLabels(inds2use);
%             fingerlabs = mainSig(vv).ResponseLabels(inds2use);

            un = unique(fingerlabs);
            num_each = sum(repmat(fingerlabs,1,numel(un))==repmat(un',numel(fingerlabs),1),1);
%             assert(all(num_each==num_each(1)))

%             runlabs = mainSig(vv).RunLabels(inds2use);
            cvlabs = mainSig(vv).SessLabels(inds2use);
            nCV = numel(unique(cvlabs));
%             cvlabs = zeros(size(runlabs));
%    
%             nCV = sum(inds2use)/2/nTrialsPerRun;
%             for mm=1:2
%                 cvlabs(maplabs==mm) = repelem(1:nCV, nTrialsPerRun);
%             end
        
            predlabs = nan(size(fingerlabs));
            
            for cv=1:nCV

                trninds = cvlabs~=cv;
                tstinds = cvlabs==cv;
                
                trnlabs = fingerlabs(trninds);
   
                un = unique(trnlabs);
                num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
%                 assert(all(num_each==num_each(1)))
                
%                  label = classify(dat_thispair(tstinds,:),dat_thispair(trninds,:),pairlabs,'diaglinear');
                
                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
    %             m = fitcsvm(dat(trninds,:), tasklabs(trninds));
    %             label = predict(m, dat(tstinds,:));
    
                predlabs(tstinds) = label;

            end
            
            assert(~any(isnan(predlabs)))

            acc = mean(predlabs==fingerlabs);
            fprintf('S%s, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},acc*100);
            all_acc(ss,vv,tt) = acc;

        end

        

    end
    
    tt=nTasks+1;
    
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    for vv=1:nVOIs
      
        inds2use = ones(size(repSig(vv).RunLabels))==1;
%   
        dat = repSig(vv).dat_avg(inds2use,:);
        dat = dat - repmat(mean(dat,2),1,size(dat,2));

        fingerlabs = nan(size(repSig(vv).IsRepeatLabels(inds2use)));
        maplabs = repSig(vv).MapLabels;
        isreplabs = repSig(vv).IsRepeatLabels(inds2use);
        
        fingerlabs(maplabs==2) = 1+isreplabs(maplabs==2);
        fingerlabs(maplabs==1) = 2-isreplabs(maplabs==1);
        
        un = unique(fingerlabs);
        num_each = sum(repmat(fingerlabs,1,numel(un))==repmat(un',numel(fingerlabs),1),1);
        assert(all(num_each==num_each(1)))

%             runlabs = mainSig(vv).RunLabels(inds2use);
        cvlabs = repSig(vv).SessLabels(inds2use);
        nCV = numel(unique(cvlabs));
%             cvlabs = zeros(size(runlabs));
%    
%             nCV = sum(inds2use)/2/nTrialsPerRun;
%             for mm=1:2
%                 cvlabs(maplabs==mm) = repelem(1:nCV, nTrialsPerRun);
%             end

        predlabs = nan(size(fingerlabs));

        for cv=1:nCV

            trninds = cvlabs~=cv;
            tstinds = cvlabs==cv;

            trnlabs = fingerlabs(trninds);

            un = unique(trnlabs);
            num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
            assert(all(num_each==num_each(1)))

%                  label = classify(dat_thispair(tstinds,:),dat_thispair(trninds,:),pairlabs,'diaglinear');

            [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
%             m = fitcsvm(dat(trninds,:), tasklabs(trninds));
%             label = predict(m, dat(tstinds,:));

            predlabs(tstinds) = label;

        end

        assert(~any(isnan(predlabs)))

        acc = mean(predlabs==fingerlabs);
        fprintf('S%s, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},acc*100);
        all_acc(ss,vv,tt) = acc;

       
    end
    
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)
%% plot decoding averaged over tasks
v2plot = [1:5,10,11,6:9,12:17];

mean_acc = squeeze(mean(mean(all_acc(:,v2plot,:),3),1))';

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],1/2,[0.1, 1],ROI_names(v2plot),[],'Accuracy',[],[])
set(gcf,'Position',[200,200,1000,400])


%% plot decoding in each task

v2plot = [1:5,10,11,6:9,12:17];

mean_acc = squeeze(mean(all_acc(:,v2plot,:),1));

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],1/2,[0.1, 1],ROI_names(v2plot),task_names,'Accuracy',[],task_colors)
set(gcf,'Position',[200,200,1000,400])
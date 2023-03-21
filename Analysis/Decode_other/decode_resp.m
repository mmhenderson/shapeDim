% decoding "response" i.e. the meaning of the response, not the finger.

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1));
% savepath = fullfile(root, 'Analysis','Decoding_results','BinaryAcc_SVM_all.mat');
savepath = fullfile(root, 'Analysis','Decoding_results','Decode_resp_all.mat');

% get ready for parallel pool operations
numcores = 8;
if isempty(gcp('nocreate'))
    parpool(numcores);
end
rndseed = 567569;

%% setup the grid

nTasks = 3;
task_names = {'Task: Linear (1)','Task: Linear (2)','Task: Checker','Repeat Detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);


nRunsTotal=nTasks*2*2*3;
nTrialsPerRun = 48;

%%
for ss=1:nSubj
   
    % First load main task data
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks+1);
        all_dprime = zeros(nSubj, nVOIs, nTasks+1);
    end

    for vv=1:nVOIs
        
        % loop over tasks, do response decoding in each separately
        for tt=1:nTasks
 
            inds2use = mainSig(vv).BoundLabels==tt;
   
            dat = mainSig(vv).dat_avg(inds2use,:);
            dat = dat - repmat(mean(dat,2),1,size(dat,2));
  
            % "catlabels" is different depending on the response mapping
            % - want to get back to "physical" categories which correspond
            % to particular parts of the shape space.
            resplabs = nan(size(mainSig(vv).CatLabels(inds2use)));
            maplabs = mainSig(vv).MapLabels(inds2use);
            catlabs = mainSig(vv).CatLabels(inds2use);
            
            resplabs(maplabs==2) = 3-catlabs(maplabs==2);
            resplabs(maplabs==1) = catlabs(maplabs==1);


            un = unique(resplabs);
            num_each = sum(repmat(resplabs,1,numel(un))==repmat(un',numel(resplabs),1),1);
            assert(all(num_each==num_each(1)))

%             cvlabs = mainSig(vv).SessLabels(inds2use);
            cvlabs = mainSig(vv).RunLabels(inds2use);
            uncv = unique(cvlabs);
            nCV = numel(uncv);

            predlabs = nan(size(resplabs));
            
            for cv=1:nCV

                trninds = cvlabs~=uncv(cv);
                tstinds = cvlabs==uncv(cv);
                
                trnlabs = resplabs(trninds);
   
                un = unique(trnlabs);
                num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
                assert(all(num_each==num_each(1)))

                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);

                predlabs(tstinds) = label;

            end
            
            assert(~any(isnan(predlabs)))

            acc = mean(predlabs==resplabs);
            dprime = get_dprime(predlabs,resplabs);
            fprintf('S%s, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},acc*100);
            all_acc(ss,vv,tt) = acc;
            all_dprime(ss,vv,tt) = dprime;

        end

        

    end
    
    tt=nTasks+1;
    
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    for vv=1:nVOIs
      
        first_trials = diff([0;repSig(vv).RunLabels])==1;
        inds2use = ~first_trials;
   
        dat = repSig(vv).dat_avg(inds2use,:);
        dat = dat - repmat(mean(dat,2),1,size(dat,2));

        resplabs = repSig(vv).IsRepeatLabels(inds2use);

        un = unique(resplabs);
        num_each = sum(repmat(resplabs,1,numel(un))==repmat(un',numel(resplabs),1),1);

%         cvlabs = repSig(vv).SessLabels(inds2use);
       cvlabs = repSig(vv).RunLabels(inds2use);
        uncv = unique(cvlabs);
        nCV = numel(uncv);

        niters = 1000;
        predlabs = nan(size(resplabs,1), niters);

        for cv=1:nCV

            trninds = find(cvlabs~=uncv(cv));
            tstinds = cvlabs==uncv(cv);

            trnlabs = resplabs(trninds);

            un = unique(trnlabs);
            num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);

            [~,smaller_group] = min(num_each);
            smaller_inds = find(trnlabs==un(smaller_group));
            larger_inds = find(trnlabs==un(3-smaller_group));
            
            predlabsthiscv = nan(sum(tstinds),niters);
            resamp_inds = nan(num_each(smaller_group), niters);
            rndseed = rndseed+1;
            rng(rndseed,'twister');
            for it = 1:niters
                resamp_inds(:,it) = datasample(larger_inds, num_each(smaller_group), 'replace', false);
            end
            parfor it =1:niters
                
                trninds_resamp = trninds([smaller_inds; resamp_inds(:,it)]);
                trnlabs_resamp = resplabs(trninds_resamp);
                num_each_resamp = sum(repmat(trnlabs_resamp,1,numel(un))==repmat(un',numel(trnlabs_resamp),1),1);
                assert(all(num_each_resamp==num_each_resamp(1)));
                
                [label,~] = normEucDistClass(dat(trninds_resamp,:),dat(tstinds,:),trnlabs_resamp);

                predlabsthiscv(:,it) = label;
            end
            predlabs(tstinds,:) = predlabsthiscv;
        end

        assert(~any(isnan(predlabs(:))))

        acc = mean(mean(predlabs==resplabs,1));
        d = nan(niters,1);
        for it = 1:niters
            d(it) = get_dprime(predlabs(:,it),resplabs);
        end
        dprime = mean(d);
        fprintf('S%s, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},acc*100);
        all_acc(ss,vv,tt) = acc;
        all_dprime(ss,vv,tt) = dprime;


       
    end
    
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)

%% plot decoding acc in each task

v2plot = [1:5,10,11,6:9,12:17];

mean_acc = squeeze(mean(all_acc(:,v2plot,:),1));

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],1/2,[0.1, 1],ROI_names(v2plot),task_names,'Accuracy',[],task_colors)
set(gcf,'Position',[200,200,1000,400])

%% plot decoding dprime in each task

v2plot = [1:5,10,11,6:9,12:17];

mean_d = squeeze(mean(all_dprime(:,v2plot,:),1));

%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
se_d = nan(size(mean_d));

plot_barsAndStars(mean_d,se_d,[],[],0,[-0.4, 1],ROI_names(v2plot),task_names,'d-prime','Binary response information (not finger)',task_colors)
set(gcf,'Position',[200,200,1000,400])
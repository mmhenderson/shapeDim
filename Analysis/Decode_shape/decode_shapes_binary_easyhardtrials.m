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
savepath = fullfile(root, 'Analysis','Decoding_results','BinaryAcc_easyhardtrials_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Linear (1)','Linear (2)','Checker'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);


nRunsTotal = nTasks*2*2*3;
nQuadrants = 4;

nDiffLevels = 2;
diff_levels = [1,0];
diff_names = {'easy','hard'};
   
quad_groups = {[1, 4; 2, 3],...
                [1, 2; 3, 4],...
                [1, 3; 2, 4]};
%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks, nDiffLevels);
    end

    for vv=1:nVOIs
        
        for tt = 1:nTasks   % loop over what task the subject was doing
            
            for dd = 1:nDiffLevels

               % difficulty means was it part of the "main" grid or not.
                difflabs = mainSig(vv).IsMainLabels;

                inds2use = mainSig(vv).BoundLabels == tt & difflabs==diff_levels(dd);

                if ss==1 && vv==1
                    fprintf('task %d %s has %d trials\n',tt,diff_names{dd},sum(inds2use));
                end
               
                dat = mainSig(vv).dat_avg(inds2use,:);

                quadlabs = mainSig(vv).QuadrantLabels(inds2use,:);

                % always decoding whichever boundary was currently active   - just
                % going to separate easy/hard.

                pairlabs = zeros(size(quadlabs));
                pairlabs(ismember(quadlabs, quad_groups{tt}(1,:)))= 1;
                pairlabs(ismember(quadlabs, quad_groups{tt}(2,:)))= 2;

                % "catlabels" is different depending on the response mapping
                % - want to get back to "physical" categories which correspond
                % to particular parts of the shape space.
                resplabs = nan(size(mainSig(vv).CatLabels(inds2use)));
                maplabs = mainSig(vv).MapLabels(inds2use);
                catlabs = mainSig(vv).CatLabels(inds2use);

                resplabs(maplabs==2) = 3-catlabs(maplabs==2);
                resplabs(maplabs==1) = catlabs(maplabs==1);

                assert(all(pairlabs==resplabs));
                
                % cross validate leaving out one run at a time
%                 cvlabs = mainSig(vv).SessLabels(inds2use,:);
                cvlabs = mainSig(vv).RunLabels(inds2use,:);
                uncv = unique(cvlabs);
                nCV = numel(uncv);

                predlabs = nan(size(dat,1),1);

                for cv=1:nCV

                    trninds = cvlabs~=uncv(cv);
                    tstinds = cvlabs==uncv(cv);
                    
                    trnlabs = pairlabs(trninds);
                    assert(numel(unique(trnlabs))==2);
                    assert(sum(trnlabs==1)==sum(trnlabs==2));

                    [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
%                     m = fitcsvm(dat(trninds,:), pairlabs(trninds));
%                     label = predict(m, dat(tstinds,:));
                    predlabs(tstinds) = label;

                end

                assert(~any(isnan(predlabs)))

                acc = mean(predlabs==pairlabs);
                fprintf('S%s, %s, task %d, %s: acc = %.2f percent\n',sublist{ss},ROI_names{vv},tt,diff_names{dd}, acc*100);
                all_acc(ss,vv,tt,dd) = acc;

            end
        end
    end
   
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)
%% plot decoding for easy and hard separately
v2plot = [1:5,10,11,6:9,12:17];

for dd=1:2
    
    mean_acc = squeeze(mean(all_acc(:,v2plot,:,dd),1));
    
%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
    se_acc = nan(size(mean_acc));
    
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.2, 1],ROI_names(v2plot),task_names,'Accuracy',diff_names{dd},task_colors(1:3,:))
    set(gcf,'Position',[200,200,1000,400])
    
end  

%% plot decoding for each task separately
v2plot = [1:5,10,11,6:9,12:17];

for tt=1:nTasks
    
    mean_acc = squeeze(mean(all_acc(:,v2plot,tt,:),1));
    
%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
    se_acc = nan(size(mean_acc));
    
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.2, 1],ROI_names(v2plot),...
        diff_names,'Accuracy',task_names{tt},[1.4*task_colors(tt,:);0.8*task_colors(tt,:)])
    set(gcf,'Position',[200,200,1000,400])
    
end  

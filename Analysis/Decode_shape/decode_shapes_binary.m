% measure overall decodability of shapes

%%
clear
close all

% sublist = {'01','02','03','04','05','06','07'};
sublist = {'01','02','03','04','05','07'};
% skipping subject 6 for now until we decide how to handle missing runs

nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end-1));
% savepath = fullfile(root, 'Analysis','Decoding_results','BinaryAcc_SVM_all.mat');
savepath = fullfile(root, 'Analysis','Decoding_results','BinaryAcc_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Linear (1)','Linear (2)','Checker','Repeat Detection'};
addpath(genpath('/usr/local/serenceslab/maggie/mFiles/Plotting tools/'))

task_colors = viridis(5);
task_colors = task_colors(1:4,:);


nRunsTotal = nTasks*2*2*3;
nQuadrants = 4;

% three different ways to do binary decoding
nBounds = 3;
bound_names = {'Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'};
quad_groups = {[1, 4; 2, 3],...
                [1, 2; 3, 4],...
                [1, 3; 2, 4]};

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks+1, nBounds);
    end

    for vv=1:nVOIs
                       
        inds2use = mainSig(vv).IsMainLabels==1;
        assert(sum(inds2use)==32*nRunsTotal);
        
        dat = mainSig(vv).dat_avg(inds2use,:);

        quadlabs = mainSig(vv).QuadrantLabels(inds2use,:);
       
        tasklabs = mainSig(vv).BoundLabels(inds2use,:);
        
        % cross validate leaving out one run at a time
%         cvlabs = mainSig(vv).SessLabels(inds2use,:);
        cvlabs = mainSig(vv).RunLabels(inds2use,:);
        nCV = numel(unique(cvlabs));
        
        predlabs = nan(size(dat,1),1);

         
        for tt = 1:nTasks   % loop over what task the subject was doing
            for bb=1:nBounds   % loop over what type of decoding to perform

                inds_this_task = tasklabs==tt;
                
                pairlabs = zeros(size(quadlabs));
                pairlabs(ismember(quadlabs, quad_groups{bb}(1,:)))= 1;
                pairlabs(ismember(quadlabs, quad_groups{bb}(2,:)))= 2;

                predlabs(inds_this_task) = nan;
               
                for cv=1:nCV

                    trninds = inds_this_task & cvlabs~=cv;
                    tstinds = inds_this_task & cvlabs==cv;
                    if sum(tstinds)==0
                        continue
                    end
                    stimlabs_trn = pairlabs(trninds);
                    assert(numel(unique(stimlabs_trn))==2);
                    assert(sum(stimlabs_trn==1)==sum(stimlabs_trn==2));

                    [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),pairlabs(trninds));
%                     m = fitcsvm(dat(trninds,:), pairlabs(trninds));
%                     label = predict(m, dat(tstinds,:));
                    predlabs(tstinds) = label;

                end

                assert(~any(isnan(predlabs(inds_this_task))))

                acc = mean(predlabs(inds_this_task)==pairlabs(inds_this_task));
                fprintf('S%s, %s, task %d, decode bound %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},tt,bb, acc*100);
                all_acc(ss,vv,tt,bb) = acc;

            end
        end
    end
    
    tt = nTasks+1;
    fn2load = fullfile(root, 'Samples', sprintf('RepeatTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    for vv=1:nVOIs
                       
        inds2use = repSig(vv).IsMainLabels==1;
       
        dat = repSig(vv).dat_avg(inds2use,:);

        points = repSig(vv).PointLabels(inds2use,:);
        center = 2.5;       
        quadlabs = zeros(size(points,1),1);
        quadlabs(points(:,1)>center & points(:,2)>center) = 1;
        quadlabs(points(:,1)<center & points(:,2)>center) = 2;
        quadlabs(points(:,1)<center & points(:,2)<center) = 3;
        quadlabs(points(:,1)>center & points(:,2)<center) = 4;

        % cross validate leaving out one run at a time
%         cvlabs = repSig(vv).SessLabels(inds2use,:);
        cvlabs = repSig(vv).RunLabels(inds2use,:);
        
        nCV = numel(unique(cvlabs));
        
        
        for bb=1:nBounds   % loop over what type of decoding to perform

            pairlabs = zeros(size(quadlabs));
            pairlabs(ismember(quadlabs, quad_groups{bb}(1,:)))= 1;
            pairlabs(ismember(quadlabs, quad_groups{bb}(2,:)))= 2;

            predlabs = nan(size(dat,1),1);

            for cv=1:nCV

                trninds = cvlabs~=cv;
                tstinds = cvlabs==cv;
                if sum(tstinds)==0
                    continue
                end
                stimlabs_trn = pairlabs(trninds);
                assert(numel(unique(stimlabs_trn))==2);
                assert(sum(stimlabs_trn==1)==sum(stimlabs_trn==2));

                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),pairlabs(trninds));
%                     m = fitcsvm(dat(trninds,:), pairlabs(trninds));
%                     label = predict(m, dat(tstinds,:));
                predlabs(tstinds) = label;

            end

            assert(~any(isnan(predlabs)))

            acc = mean(predlabs==pairlabs);
            fprintf('S%s, %s, repeat task, decode bound %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},bb, acc*100);
            all_acc(ss,vv,tt,bb) = acc;

        end
      
    end
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)
%% plot decoding along each axis separately
v2plot = [1:5,10,11,6:9];

for bb=1:3
    
    mean_acc = squeeze(mean(all_acc(:,v2plot,:,bb),1));
    
    se_acc = squeeze(std(all_acc(:,v2plot,:,bb),[],1)./sqrt(nSubj-1));
%     se_acc = nan(size(mean_acc));
    
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names,'Accuracy',sprintf('%s',bound_names{bb}),task_colors)
    set(gcf,'Position',[200,200,1000,400])
    
end   

%% plot decoding for just the two linear tasks
v2plot = [1:5,10,11,6:9];
task2plot = 1:2;
for bb=1:2
    
    mean_acc = squeeze(mean(all_acc(:,v2plot,task2plot,bb),1));
    
    se_acc = squeeze(std(all_acc(:,v2plot,task2plot,bb),[],1)./sqrt(nSubj-1));
%     se_acc = nan(size(mean_acc));
    
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names(task2plot),'Accuracy',sprintf('%s',bound_names{bb}),task_colors(task2plot,:))
    set(gcf,'Position',[200,200,1000,400])
    
end   

% %% plot decoding along each axis separately
% v2plot = [1:5,10,11,6:9,12:17];
% 
% for tt = 1:3
%     
%     mean_acc = squeeze(mean(all_acc(:,v2plot,tt,:),1));
%     
% %     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
%     se_acc = nan(size(mean_acc));
%     
%     plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),bound_names,'Accuracy',sprintf('%s',task_names{tt}),[])
% 
% end    
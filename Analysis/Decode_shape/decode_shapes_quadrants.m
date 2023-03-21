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
savepath = fullfile(root, 'Analysis','Decoding_results','Decode_quadrant_all.mat');

%% setup the grid

nTasks = 3;
task_names = {'Task: Linear (1)','Task: Linear (2)','Task: Checker', 'Repeat detection'};
task_colors = viridis(5);
task_colors = task_colors(1:4,:);

nQuadrants = 4;
quad_pairs = combnk(1:nQuadrants,2);
nPairs = size(quad_pairs,1);

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('MainTaskSignalByTrial_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nTasks+1, nPairs);
    end

    for vv=1:nVOIs
        
        for tt = 1:nTasks   % loop over what task the subject was doing
                       
            inds2use = mainSig(vv).BoundLabels==tt & mainSig(vv).IsMainLabels==1;
        
            dat = mainSig(vv).dat_avg(inds2use,:);

            quadlabs = mainSig(vv).QuadrantLabels(inds2use,:);

            % cross validate leaving out one sess at a time
            cvlabs = mainSig(vv).SessLabels(inds2use);
            nCV = numel(unique(cvlabs));

            predlabs = nan(size(quadlabs));
            
            for pp=1:nPairs   % loop over what type of decoding to perform

                inds2use = quadlabs==quad_pairs(pp,1) | quadlabs==quad_pairs(pp,2);
                predlabs(inds2use) = nan;
                
                for cv=1:nCV

                    trninds = inds2use & cvlabs~=cv;
                    tstinds = inds2use & cvlabs==cv;
                    
                    trnlabs = quadlabs(trninds);
                    un = unique(trnlabs);
                    num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
                    assert(all(num_each==num_each(1)))

                    [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
%                     m = fitcsvm(dat(trninds,:), pairlabs(trninds));
%                     label = predict(m, dat(tstinds,:));
                    predlabs(tstinds) = label;

                end

                assert(~any(isnan(predlabs(inds2use))))

                acc = mean(predlabs(inds2use)==quadlabs(inds2use));
                fprintf('S%s, %s, task %d, pair %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},tt,pp, acc*100);
                all_acc(ss,vv,tt,pp) = acc;

            end
        end
    end
    
    tt=nTasks+1;
    
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

        % cross validate leaving out one sess at a time
        cvlabs = repSig(vv).SessLabels(inds2use);
        nCV = numel(unique(cvlabs));

        predlabs = nan(size(quadlabs));

        for pp=1:nPairs   % loop over what type of decoding to perform

            inds2use = quadlabs==quad_pairs(pp,1) | quadlabs==quad_pairs(pp,2);
            predlabs(inds2use) = nan;

            for cv=1:nCV

                trninds = inds2use & cvlabs~=cv;
                tstinds = inds2use & cvlabs==cv;

                trnlabs = quadlabs(trninds);
                un = unique(trnlabs);
                num_each = sum(repmat(trnlabs,1,numel(un))==repmat(un',numel(trnlabs),1),1);
                assert(all(num_each==num_each(1)))

                [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),trnlabs);
%                     m = fitcsvm(dat(trninds,:), pairlabs(trninds));
%                     label = predict(m, dat(tstinds,:));
                predlabs(tstinds) = label;

            end

            assert(~any(isnan(predlabs(inds2use))))

            acc = mean(predlabs(inds2use)==quadlabs(inds2use));
            fprintf('S%s, %s, repeat task, pair %d: acc = %.2f percent\n',sublist{ss},ROI_names{vv},pp, acc*100);
            all_acc(ss,vv,tt,pp) = acc;

        end
      
    end
    
end

%%
save(savepath, 'all_acc','ROI_names');

%%
load(savepath)
%% plot decoding along each axis separately
v2plot = [1:5,10,11,6:9,12:17];

for tt=1:4
    
    mean_acc = squeeze(mean(mean(all_acc(:,v2plot,tt,:),4),1))';
    
%     se_acc = squeeze(std(all_acc(:,v2plot,tt,:),[],1)./sqrt(nSubj-1));
    se_acc = nan(size(mean_acc));
    
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),[],'Accuracy',sprintf('Decode quadrant, %s',task_names{tt}),task_colors(tt,:))
    set(gcf,'Position',[200,200,1000,400])
end   

%% plot decoding along each axis separately
v2plot = [1:5,10,11,6:9,12:17];

mean_acc = squeeze(mean(mean(all_acc(:,v2plot,:,:),4),1));

se_acc = nan(size(mean_acc));

plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],ROI_names(v2plot),task_names,'Accuracy','Decode quadrant',task_colors)

set(gcf,'Position',[200,200,1000,400])
  
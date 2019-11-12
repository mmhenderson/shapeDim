% measure overall decodability of shapes
% downsampling trials here to see if we can still do ok with fewer trials

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end)-1);
savepath = fullfile(root, 'Analysis','Decoding_results','PairwiseAcc_Retinotopic_downsample.mat');

%% setup the grid
nSets = 2;
nStimsPerSet = 16;
pairs = combnk(1:nStimsPerSet, 2);
nPairs = size(pairs,1);
     
pairwise_dist = zeros(size(pairs,1),1);

[xpoints,ypoints] = meshgrid(round(linspace(0.2, 4.8, 4),1), round(linspace(0.2, 4.8, 4),1));
xpoints = xpoints(:);
ypoints = ypoints(:);
points = [xpoints, ypoints];

which_dim = zeros(size(points,1),1);
for pp=1:size(pairs,1)
    pairwise_dist(pp) = sqrt(sum((points(pairs(pp,1),:) - points(pairs(pp,2),:)).^2));
    
    % do the points share the same y coordinate? if so, the discrimination
    % between them is along the x-coordinate.
    if points(pairs(pp,1),2)==points(pairs(pp,2),2)
        which_dim(pp) = 1;
    elseif points(pairs(pp,1),1)==points(pairs(pp,2),1)
        which_dim(pp) = 2;
    end
end

pairwise_dist_round = round(pairwise_dist*2)./2;

[unvals, ia,ib] = unique(pairwise_dist_round);

numeach = sum(repmat(ib, 1, numel(unvals))==repmat((1:numel(unvals)), numel(ib,1), 1), 1);

%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('NeutralTaskSignalByTrial_Retinotopic_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        all_acc = zeros(nSubj, nVOIs, nSets, nPairs);
    end

    for vv=1:nVOIs
        
        dat = mainSig(vv).mainDat;
        
        runlabs = mainSig(vv).runLabs;
        nRuns = numel(unique(runlabs));
        trials2use = runlabs<=nRuns/2;
        nTrialsPerRun = numel(runlabs)/nRuns;
        
        stimlabs = mainSig(vv).StimLabels;
        
        % taking just first half of teh data
        dat = dat(trials2use,:);
        runlabs=  runlabs(trials2use,:);
        stimlabs = stimlabs(trials2use,:);
        
        nRuns = nRuns/2;
        
        % sets are always alternating on even and odd runs
        setlabs = mod(runlabs,2)+1;
        
        % cross validate leaving out one run at a time
        cvlabs = repelem(1:nRuns/2, nTrialsPerRun*2)';
        nCV = numel(unique(cvlabs));
        
        for st = 1:nSets
            for pp=1:nPairs

                pairinds = (stimlabs==pairs(pp,1) | stimlabs==pairs(pp,2)) & setlabs==st;

                label_all = nan(size(dat,1),1);

                for cv=1:nCV

                    trninds = pairinds & cvlabs~=cv;
                    tstinds = pairinds & cvlabs==cv;

                    stimlabs_trn = stimlabs(trninds);
                    assert(numel(unique(stimlabs_trn))==2);
                    assert(sum(stimlabs_trn==pairs(pp,1))==sum(stimlabs_trn==pairs(pp,2)));

                    [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),stimlabs(trninds));

                    label_all(tstinds) = label;

                end

                assert(~any(isnan(label_all(pairinds))))

                acc = mean(label_all(pairinds)==stimlabs(pairinds));
                fprintf('S%s, %s, set %d, pair %d: acc = %.2f percent\n',sublist{ss},my_areas{vv},st,pp, acc*100);
                all_acc(ss,vv,st,pp) = acc;

            end
        end
    end
end

%%
save(savepath, 'all_acc','my_areas');

%%
for st=1:nSets
    mean_acc =  mean(squeeze(all_acc(:,:,st,:)), 2);
    se_acc = std(squeeze(all_acc(:,:,st,:)),[],2)./sqrt(nPairs);

    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],my_areas,[],'Accuracy',sprintf('Set %d\nRetinotopic ROIs: average of all pairwise classifiers',st),[])
end
%%

v2plot = [1:6];
ss=1;

cols = plasma(numel(v2plot)+1);
cols = cols(1:numel(v2plot),:);

for st=1:nSets

    figure;hold all;

    for vv=1:numel(v2plot)


        meanvals = zeros(size(unvals));
        sevals = zeros(size(unvals));

        for uu=1:length(unvals)

            inds = pairwise_dist_round==unvals(uu);
            meanvals(uu) = mean(all_acc(ss,v2plot(vv),st,inds));
            sevals(uu) = std(all_acc(ss,v2plot(vv),st,inds))./sqrt(numel(inds));

        end

        errorbar(unvals, meanvals, sevals,'Color',cols(vv,:),'LineWidth',2);

        xlabel('Distance between shapes in feature space')
        ylabel('Decoding accuracy');


    end
    line(get(gca,'XLim'),[0.5, 0.5],'Color','k')

    title(sprintf('Set %d\nDecoding accuracy vs. distance in feature space',st));
    legend(my_areas(v2plot),'Location','EastOutside');
    set(gcf,'Color','w');

    ylim([0.2, 1])

end

%% plot decoding along each axis separately
v2plot = [1:6];
for st = 1:2
    
    mean_vals = zeros(nSubj, numel(v2plot), 3);
    
    for dd=1:3

        pairs2use = which_dim+1==dd;
        mean_vals(:,:,dd) = squeeze(mean(all_acc(:,v2plot,st,pairs2use),4));
        se_vals(:,:,dd) = squeeze(std(all_acc(:,v2plot,st,pairs2use),[],4))./sqrt(sum(pairs2use));
        
    end

%     if nSubj==1
        mean_acc = squeeze(mean_vals);
        se_acc = squeeze(se_vals);
%     else
%         mean_acc =  mean(mean_vals,1);
%         se_acc = std(mean_vals,[],1)./sqrt(nSubj);
%     end
    plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],my_areas,{'mixed','dim 1','dim 2'},'Accuracy',sprintf('Set %d\nRetinotopic ROIs: average of all pairwise classifiers',st),[])

end    
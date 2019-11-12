% measure overall decodability of shapes

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end)-1);
savepath = fullfile(root, 'Analysis','Decoding_results','PairwiseAcc_Retinotopic_collapseDims.mat');

%% setup the grid
nSets = 2;   
nAxes = 2;
nStimsPerAxis = 4;

pairs = combnk(1:nStimsPerAxis, 2);
nPairs = size(pairs,1);
 
pairwise_dist = zeros(size(pairs,1),1);

points = round(linspace(0.2, 4.8, 4),1);

for pp=1:size(pairs,1)
    pairwise_dist(pp) = abs(points(pairs(pp,1)) - points(pairs(pp,2)));
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
        all_acc = zeros(nSubj, nVOIs, nSets, nAxes, nPairs);
    end

    for vv=1:nVOIs
        
        dat = mainSig(vv).mainDat;
        
        runlabs = mainSig(vv).runLabs;
        nRuns = numel(unique(runlabs));
        nTrialsPerRun = numel(runlabs)/nRuns;
        
        
        % sets are always alternating on even and odd runs
        setlabs = mod(runlabs,2)+1;
        
        % cross validate leaving out one run at a time
        cvlabs = repelem(1:nRuns/2, nTrialsPerRun*2)';
        nCV = numel(unique(cvlabs));
        
        for st = 1:nSets
            
            for ax = 1:nAxes
                
                % take out just the labels for my axis of interest,
                % ignoring values of the other feature
                % stim labels are 1-4 integers 
                stimlabs = mainSig(vv).ShapeCoords(:,ax);
                [u,ia,ib] = unique(stimlabs);
                assert(all(round(u',1)==points));
                stimlabs = ib;
                
                for pp=1:nPairs

                    pairinds = (stimlabs==pairs(pp,1) | stimlabs==pairs(pp,2)) & setlabs==st;

                    label_all = nan(size(dat,1),1);

                    for cv=1:nCV

                        trninds = pairinds & cvlabs~=cv;
                        tstinds = pairinds & cvlabs==cv;

                        stimlabs_trn = stimlabs(trninds);
                        assert(numel(unique(stimlabs_trn))==2);
                        % double check the set is balanced for two
                        % categories
                        assert(sum(stimlabs_trn==pairs(pp,1))==sum(stimlabs_trn==pairs(pp,2)));

                        [label,~] = normEucDistClass(dat(trninds,:),dat(tstinds,:),stimlabs(trninds));

                        label_all(tstinds) = label;

                    end

                    assert(~any(isnan(label_all(pairinds))))

                    acc = mean(label_all(pairinds)==stimlabs(pairinds));
                    fprintf('S%s, %s, set %d, axis %d, pair %d: acc = %.2f percent\n',sublist{ss},my_areas{vv},st,ax,pp,acc*100);
                    all_acc(ss,vv,st,ax,pp) = acc;

                end
            end
        end
    end
end

%%
save(savepath, 'all_acc','my_areas');

%%
for st=1:nSets
    for ax = 1:nAxes
        mean_acc =  mean(squeeze(all_acc(:,:,st,ax,:)), 2);
        se_acc = std(squeeze(all_acc(:,:,st,ax,:)),[],2)./sqrt(nPairs);

        plot_barsAndStars(mean_acc,se_acc,[],[],0.5,[0.4, 1],my_areas,[],'Accuracy',sprintf('Set %d, Axis %d\nRetinotopic ROIs: average of all pairwise classifiers',st,ax),[])
    end
end
%%

v2plot = [1:6];
ss=1;

cols = plasma(numel(v2plot)+1);
cols = cols(1:numel(v2plot),:);

for st=1:nSets
    
    for ax = 1:nAxes

        figure;hold all;

        for vv=1:numel(v2plot)


            meanvals = zeros(size(unvals));
            sevals = zeros(size(unvals));

            for uu=1:length(unvals)

                inds = pairwise_dist_round==unvals(uu);
                meanvals(uu) = mean(all_acc(ss,v2plot(vv),st,ax,inds));
                sevals(uu) = std(all_acc(ss,v2plot(vv),st,ax,inds))./sqrt(numel(inds));

            end

            errorbar(unvals, meanvals, sevals,'Color',cols(vv,:),'LineWidth',2);

            xlabel(sprintf('Distance between shapes along axis %d',ax))
            ylabel('Decoding accuracy');


        end
        line(get(gca,'XLim'),[0.5, 0.5],'Color','k')

        title(sprintf('Set %d, Axis %d\nDecoding accuracy collapsed over other axis',st,ax));
        legend(my_areas(v2plot),'Location','EastOutside');
        set(gcf,'Color','w');

        ylim([0.2, 1])
    end
end

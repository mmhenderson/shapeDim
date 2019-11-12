% measure overall decodability of shapes

%%
clear
close all

sublist = {'01'};
nSubj = length(sublist);
my_dir = pwd;
filesepinds = find(my_dir==filesep);
root = my_dir(1:filesepinds(end)-1);
savepath = fullfile(root, 'Analysis','Decoding_results','PairwiseAcc_Retinotopic.mat');

nStims = 16;
[xpoints,ypoints] = meshgrid(round(linspace(0.2, 4.8, 4),1), round(linspace(0.2, 4.8, 4),1));
xpoints = xpoints(:);
ypoints = ypoints(:);
points = [xpoints, ypoints];
%%
for ss=1:nSubj
   
    fn2load = fullfile(root, 'Samples', sprintf('NeutralTaskSignalByTrial_Retinotopic_S%s.mat', sublist{ss}));
    load(fn2load);
    
    nVOIs = numel(mainSig);
    
    if ss==1
        euc_dist = zeros(nSubj, nVOIs, nStims, nStims);
    end

    for vv=1:nVOIs
        
        dat = mainSig(vv).mainDat;
        stimlabs = mainSig(vv).StimLabels;
        
        for ss1 = 1:nStims
            
            for ss2 = ss1+1:nStims
                
                dat1 = dat(stimlabs==ss1,:);
                dat2 = dat(stimlabs==ss2,:);
                
                dist = get_normEucDist(dat1,dat2);
                
                euc_dist(ss,vv,ss1,ss2) = dist;
                euc_dist(ss,vv,ss2,ss1) = dist;
            end
        end      
    end
end

%%
ss=1;
cols = plasma(nStims);
v2plot = 1:4;
for vv=1:numel(v2plot)
   
    D = squeeze(euc_dist(ss,v2plot(vv),:,:));
    
    reduced_rep = mdscale(D,2);
    
    figure;hold all;
    legendlabs = {};
    for pp=1:size(reduced_rep,1)
        
        plot(reduced_rep(pp,1),reduced_rep(pp,2),'o','MarkerFaceColor',cols(pp,:))
        legendlabs{pp} = sprintf('%.1f, %.1f',points(pp,1),points(pp,2));
    end
    legend(legendlabs,'Location','EastOutside')
    title(sprintf('MDS: %s\n',my_areas{v2plot(vv)}));
    xlabel('MDS axis 1');
    ylabel('MDS axis 2');
    
end
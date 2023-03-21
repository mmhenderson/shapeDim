%% find FDR threshold for all my subjects (and save out this image)

sbj = {'01'};
exp_path = '/mnt/neurosphere/serenceslab2/maggie/OriSpin/';

for n = 1:numel(sbj)
    
   
    
    % find stats map
    stats_path = [exp_path 'AnalyzeLocalizer/S' char(sbj(n)) '/feats/AllSessionsFE.gfeat/cope1.feat/stats/'];
    cd(stats_path)

     % get FE dof
    [~, my_FE_dof] = unix('fslstats tdof_t1 -M'); % gives the dof I need for fixed effects (which is the test I ran), output is a string
    
    % make log p-stats map
    unix(['ttologp -logpout logp1 varcope1 cope1 ', num2str(str2num(my_FE_dof))]);
    
    % convert from log to p-value map
    unix(['fslmaths logp1 -exp p1']);
    
    % do FDR on p-value map and get probability threshold
    [~,prob_thresh] = unix('fdr -i p1 -q 0.05');

    % go from my p-threshold back into my t-stat
    my_t(n) = abs(tinv(str2num(prob_thresh(28:end)),str2num(my_FE_dof)));
    
end

fns = [exp_path, 'Samples/t_cutoff.mat'];
save(fns,'my_t');



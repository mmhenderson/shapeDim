function dist = get_normEucDist(dat1,dat2)
% """Calculate the normalized euclidean distance (d') between the means of
%     two clouds of data points.
% 
%     Args:
%       dat1: [nPts1 x nWeights] (voxels, spike rates, neural network weights)
%       dat2: [nPts2 x nWeights] 
%      
%     Returns:
%        normEucDist (single value)
%     """

    assert(size(dat1,2)==size(dat2,2))
    
    npts1 = size(dat1,1);
    npts2 = size(dat2,1);

    var1 = var(dat1,[],1);
    var2 = var(dat2,[],1);
    
    pooled_var = (var1.*(npts1-1)+var2.*(npts2-1))./(npts1-1+npts2-1);
    
    mean1 = mean(dat1,1);
    mean2 = mean(dat2,1);
    
    sq = (mean1-mean2).^2./pooled_var;
    dist = sqrt(sum(sq));

end

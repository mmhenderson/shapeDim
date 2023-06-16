import numpy as np
import scipy.stats
import warnings


def lin_reg(x,y):
   
    if len(x.shape)==1:
        x_mat = x[:,np.newaxis]
    else:
        x_mat = x
    if len(y.shape)==1:
        y_mat = y[:,np.newaxis]
    else:
        y_mat = y
        
    n_pts = x_mat.shape[0]
    assert(y_mat.shape[0]==n_pts)
    
    X = np.concatenate([x_mat, np.ones((n_pts,1))], axis=1)
    reg_coeffs = np.linalg.pinv(X) @ y_mat
    yhat = X @ reg_coeffs
    
    actual = np.squeeze(y_mat)
    pred = np.squeeze(yhat)
    ssres = np.sum(np.power((actual - pred),2));
    sstot = np.sum(np.power((actual - np.mean(actual)),2));
    r2 = 1-(ssres/sstot)
    
    return yhat, reg_coeffs, r2

def get_dprime(predlabs,reallabs,un=None):
    """ 
    Calculate d' for predicted and actual values. Works for multiple classes.
    """

    predlabs==np.squeeze(predlabs)
    reallabs==np.squeeze(reallabs)
    if len(predlabs)!=len(reallabs):
        raise ValueError('real and predicted labels do not match')
    if len(predlabs.shape)>1 or len(reallabs.shape)>1:
        raise ValueError('need to have 1d inputs')
    if un is None:
        un = np.unique(reallabs)
    if not np.all(np.isin(np.unique(predlabs), un)):
        print('Warning: some labels in pred are not included in real labels! Will return nan')
        return np.nan
    
    hrz=np.zeros((len(un),1));
    fpz=np.zeros((len(un),1));

    n_trials = len(predlabs);

    #loop over class labels, get a hit rate and false pos for each (treating
    #any other category as non-hit)
    for ii in range(len(un)):

        if np.sum(reallabs==un[ii])==0 or np.sum(reallabs!=un[ii])==0:

            # if one of the categories is completely absent - this will return a
            # nan dprime value
            return np.nan

        else:

            hr = np.sum((predlabs==un[ii]) & (reallabs==un[ii]))/np.sum(reallabs==un[ii]);
            fp = np.sum((predlabs==un[ii]) & (reallabs!=un[ii]))/np.sum(reallabs!=un[ii]);    

            # make sure this never ends up infinite
            # correction from Macmillan & Creelman, use 1-1/2N or 1/2N in place
            # of 1 or 0 
            if hr==0:
                hr=1/(2*n_trials)
            if fp==0:
                fp=1/(2*n_trials)
            if hr==1:
                hr=1-1/(2*n_trials)
            if fp==1:
                fp=1-1/(2*n_trials);

        # convert to z score (this is like percentile - so 50% hr would be zscore=0)
        hrz[ii]=scipy.stats.norm.ppf(hr,0,1);
        fpz[ii]=scipy.stats.norm.ppf(fp,0,1);

    # dprime is the mean of individual dprimes (for two classes, they will be
    # same value)
    dprime = np.mean(hrz-fpz);

    return dprime




def compute_partial_corr(x, y, c, return_p=False):

    """
    Compute the partial correlation coefficient between x and y, 
    controlling for the variables in covariates "c". 
    Uses linear regression based method.
    Inputs: 
        x [n_samples,] or [n_samples,1]
        y [n_samples,] or [n_samples,1]
        c [n_samples,] or [n_samples,n_covariates]
        
    Outputs:
        partial_corr, a single value for the partial correlation coefficient.
    """
    
    if len(x.shape)==1:
        x = x[:,np.newaxis]        
    if len(y.shape)==1:
        y = y[:,np.newaxis]
    if len(c.shape)==1:
        c = c[:,np.newaxis]
    n_trials = x.shape[0]
    assert(y.shape[0]==n_trials and c.shape[0]==n_trials)
    
    # first predict x from the other vars
    model1_preds = np.concatenate([c, np.ones((n_trials,1))], axis=1)
    model1_coeffs = np.linalg.pinv(model1_preds) @ x
    model1_yhat = model1_preds @ model1_coeffs
    model1_resids = model1_yhat - x
   
    # then predict y from the other vars
    model2_preds = np.concatenate([c, np.ones((n_trials,1))], axis=1)
    model2_coeffs = np.linalg.pinv(model2_preds) @ y
    model2_yhat = model2_preds @ model2_coeffs
    model2_resids = model2_yhat - y

    # correlate the residuals to get partial correlation.
    if return_p:
        partial_corr, p = scipy.stats.pearsonr(model1_resids[:,0], model2_resids[:,0])
        return partial_corr, p
    else:
        partial_corr = numpy_corrcoef_warn(model1_resids[:,0], model2_resids[:,0])[0,1]
        return partial_corr
   
    

def compute_partial_corr_formula(x,y,c):
 
    """
    Code to compute the partial correlation between x and y, controlling
    for covariate c. 
    Based on the correlation coefficients between each pair of variables.
    Also computes estimated standardized beta weight. 
    """
    x = np.squeeze(x); 
    y = np.squeeze(y);
    c = np.squeeze(c);
    ryx = np.corrcoef(x,y)[0,1]
    ryc = np.corrcoef(c,y)[0,1]
    rxc = np.corrcoef(x,c)[0,1]

    # partial correlation coefficient
    partial_corr = (ryx - ryc*rxc)/np.sqrt((1-ryc**2)*(1-rxc**2))
  
    # equivalent to standardized beta weight from a multiple linear regression
    # would be set up like [x, c, intercept] @ w = y
    # this is the weight for x.
    # standarized beta = raw beta * std(x)/std(y)
    beta = (ryx - ryc*rxc)/(1-rxc**2)
    
    return partial_corr, beta

# Some functions that wrap basic numpy/scipy functions, but will print 
# more useful warnings when a problem arises

def numpy_corrcoef_warn(a,b):
    
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            cc = np.corrcoef(a,b)
        except RuntimeWarning as e:
            print('Warning: problem computing correlation coefficient')
            print('shape a: ',a.shape)
            print('shape b: ',b.shape)
            print('sum a: %.9f'%np.sum(a))
            print('sum b: %.9f'%np.sum(b))
            print('std a: %.9f'%np.std(a))
            print('std b: %.9f'%np.std(b))
            print(e)
            warnings.filterwarnings('ignore')
            cc = np.corrcoef(a,b)
            
    if np.any(np.isnan(cc)):
        print('There are nans in correlation coefficient')
    
    return cc



def paired_ttest_nonpar(vals1, vals2, n_iter=1000, rndseed=None):
    
    if rndseed is None:
        rndseed = int(time.strftime('%M%H%d', time.localtime()))
    np.random.seed(rndseed)
        
    real_diff = np.mean(vals1-vals2)    
    
    shuff_diffs = np.zeros((n_iter,))
    
    for ii in range(n_iter):
        
        shuff_vals = np.array([vals1, vals2])
        # randomly swap the positions of values within a pair, with 50% prob
        which_swap = np.random.normal(0,1,[len(vals1),])>0
        shuff_vals[:,which_swap] = np.flipud(shuff_vals[:, which_swap])
    
        shuff_diffs[ii] = np.mean(shuff_vals[0,:]-shuff_vals[1,:])
    
    # pvalue for two-tailed test
    pval_twotailed = np.minimum( np.mean(shuff_diffs<=real_diff), \
                                 np.mean(shuff_diffs>=real_diff)) * 2
    
    return pval_twotailed, real_diff

import numpy as np
import os, sys
import pandas as pd
import scipy.stats
import scipy.io as spio

from code_utils import file_utils

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/


def get_hrfs_maintask(ss, nTRs_to_model = 20, use_bigIPS = True, concat_IPS = True):
    
    if use_bigIPS:
        sample_fn = os.path.join(root, 'Samples','SampleFile_bigIPS_S%02d.mat'%ss)
    else:
        sample_fn = os.path.join(root, 'Samples','SampleFile_S%02d.mat'%ss)


    sample_fn = os.path.join(root, 'Samples','SampleFile_bigIPS_S%02d.mat'%ss)

    samples = file_utils.load_samplefile_h5py(sample_fn)

    samples['all_vox_concat'].shape
    
    # load the "timing" file (what happened on each TR)
    # this is made in GetEventTiming.m
    timing_fn = os.path.join(root, 'Samples','TimingFile_S%02d.mat'%ss)
    # print('loading from %s'%timing_fn)
    timing = spio.loadmat(timing_fn, squeeze_me=True, struct_as_record=False)
    main = file_utils._todict(timing['main'])
    rep = file_utils._todict(timing['rep'])
    
    # load the behav data structure, this is [n_trials,] and has more trial attributes
    behav_fn = os.path.join(root, 'DataBehavior', 'S%02d'%ss, \
                                          'S%02d_maintask_preproc_all.csv'%ss)
    # print('loading from %s'%behav_fn)
    bdat = pd.read_csv(behav_fn, index_col=0)
    
    roi_names = samples['ROI_names']
    n_rois = len(roi_names)
    n_hemis = len(samples['hemis'])
    
    h = []
    
    for rr in range(n_rois):

        dat_this_roi = []
        
        for hh in range(n_hemis):

            roi_inds_num = np.array(samples['ROIs']['voxel_inds'][rr][hh])

            inds_this_hemi = np.isin(np.array(samples['all_vox_concat']), roi_inds_num)[:,0]

            if np.sum(inds_this_hemi)>0:
                # print(np.sum(inds_this_hemi))
                # samples['samplesMain'] is originally [n_vox x nTRs]
                # transpose because i want [nTRs x n_vox]
                dat_this_hemi = samples['samplesMain'][inds_this_hemi,:].T
                dat_this_roi += [dat_this_hemi]
                
        dat_this_roi = np.concatenate(dat_this_roi, axis=1)
        
        # count things, check the counts
        nTRs = 327-16;
        nTRs_total, n_vox = dat_this_roi.shape
        assert(np.mod(nTRs_total, nTRs)==0)
        n_runs = int(nTRs_total / nTRs)
        
        run_labels = np.array(main['RunLabels'])
        assert(n_runs==len(np.unique(run_labels)))
        assert(nTRs_total==len(run_labels))
        
        # zscore the data from each run to normalize.
        for run in np.unique(run_labels):

            run_inds = run_labels==run
            assert(np.sum(run_inds)==nTRs)

            dat_this_roi[run_inds,:] = scipy.stats.zscore(dat_this_roi[run_inds,:], axis=0)
        
        
        # label the onset of each trial
        # 1 = stim on, 0 = stim off
        event_diff  = np.diff(np.array([0] + main['EventLabels']))
        event_diff_reshaped = np.reshape(event_diff, [nTRs, n_runs], order='F')
        trial_onset_bool = event_diff_reshaped==1;
        trial_onset_bool = np.reshape(trial_onset_bool, [nTRs*n_runs,1], order='F')
        trial_onset_num = np.where(trial_onset_bool)[0]

        n_trials_per_run = 48
        n_trials = n_runs*n_trials_per_run
        assert(len(trial_onset_num)==n_trials)
        assert(np.all(np.unique(bdat['run_overall'])==np.unique(run_labels)))

        # defining the trial conditions we want HRF for
        task_names = ['Linear-1','Linear-2','Checker']
        trial_type_names = ['Easy','Hard']

        n_tasks = 3;
        n_trial_types = 2;
        # n_trial_conds = n_tasks * n_trial_types + 1
        n_trial_conds = n_tasks * n_trial_types

        is_main_grid = np.array(bdat['is_main_grid']==1)
        task_labs = np.array(bdat['task'])

        trial_cond_names = []

        # making a matrix of onsets for each trial cond
        trial_cond_onsets = np.zeros((nTRs_total,n_trial_conds))

        ci = -1
        for ti,tt in enumerate([1,2,3]):

            for tyi in range(n_trial_types):

                ci+=1

                # figure out which of the onsets belong to this trial condition
                if tyi==0:
                    trial_inds = (task_labs==tt) & is_main_grid
                else:
                    trial_inds = (task_labs==tt) & ~is_main_grid

                onsets = trial_onset_num[trial_inds].astype(int)

                trial_cond_onsets[onsets,ci] = 1

                trial_cond_names += ['%s %s'%(task_names[ti], trial_type_names[tyi])]

#         # add response as another predictor here
#         ci+=1
#         resp_times = np.array(bdat['rt'])
#         resp_times_trs = np.round(resp_times/0.8)

#         resp_onsets = trial_onset_num + resp_times_trs
#         resp_onsets = resp_onsets[~np.isnan(resp_onsets)].astype(int)

#         trial_cond_onsets[resp_onsets,ci] = 1
#         trial_cond_names += ['Response']

        # do the deconvolution here
        
        hrfs, b, vox_r2 = do_deconv(dat_this_roi, trial_cond_onsets, run_labels, nTRs_to_model=nTRs_to_model)
        
        hroi = dict()
        
        hroi['hrfs'] = hrfs
        hroi['b'] = b
        hroi['vox_r2'] = vox_r2
        
        h += [hroi]
        
    if concat_IPS:

        # going to combine all the 4 ips subregions together
        # attempting to boost signal since these areas are small in some subjects.
        # the rois will now go:
        # ['V1','V2','V3','V3AB','hV4','LO1','LO2','IPSall']
        ips_roi = [5,6,7,8]

        roi_names = roi_names[0:5] + roi_names[9:] + ['IPSall']
        n_rois = len(roi_names)

        h_new = h[0:5] + h[9:]
        
        newips = dict()
        newips['hrfs'] = np.concatenate([h[rr]['hrfs'] for rr in ips_roi], axis=2)
        newips['b'] = np.concatenate([h[rr]['b'] for rr in ips_roi], axis=1)
        newips['vox_r2'] = np.concatenate([h[rr]['vox_r2'] for rr in ips_roi], axis=0)
        
        h_new += [newips]
        
        h = h_new
            
    return h, trial_cond_names, roi_names



def do_deconv(data, trial_cond_onsets, run_labels, nTRs_to_model=14):
    
    """do deconvolution for fmri data
    loosely based on a matlab function doDecon_forRR.m 
    that i copied from JS/RR"""

    nTRs_total, n_vox = data.shape
    n_trial_conds = trial_cond_onsets.shape[1]
    assert(trial_cond_onsets.shape[0]==nTRs_total)
    assert(len(run_labels)==nTRs_total)
    n_runs = len(np.unique(run_labels))
    nTRs = nTRs_total/n_runs
    
    # there are predictors for each TR of each condition, and for each run
    n_conds_total = n_trial_conds * nTRs_to_model + n_runs 

    X = np.zeros((nTRs_total, n_conds_total+1))

    for ri, rr in enumerate(np.unique(run_labels)):

        run_inds = (run_labels==rr)
        assert(np.sum(run_inds)==nTRs)

        # make predictors for each condition
        # these go in same columns as with other runs
        pred_ind = 0
        for ci in range(n_trial_conds):

            # onsets of all events this condition, this run
            cond_onsets = np.where((trial_cond_onsets[:,ci]==1) & run_inds)[0]

            for tr in range(nTRs_to_model):

                # shift by one TR
                onsets_adj = cond_onsets + tr

                # make sure we are still in current run, not bleeding into next
                onsets_adj = onsets_adj[onsets_adj<nTRs_total]
                onsets_adj = onsets_adj[run_inds[onsets_adj]]
                # assert(np.all(run_inds[onsets_adj]))

                # print(onsets_adj)
                X[onsets_adj,pred_ind] = 1;

                pred_ind+=1

        # make constant predictor for this run
        pred_ind = n_trial_conds * nTRs_to_model + ri
        X[run_inds, pred_ind] = 1
        
    # add a column of all 1s, intercept over all runs
    X[:, n_conds_total] = 1;

    # check that X is full-rank, not degenerate
    # except for the last column which is constant ones
    assert(np.linalg.matrix_rank(X) >= (X.shape[1]-1))

    # solve for beta weights
    # b is [n_conds_total x n_vox]
    b = np.linalg.pinv(X) @ data

    # reshape to [nTRs x n_trial_conds x n_vox]
    # this removes the constant-run predictors
    hrfs = np.reshape(b[0:n_trial_conds * nTRs_to_model, :], [nTRs_to_model, n_trial_conds, n_vox], order='F')

    run_const = b[n_trial_conds * nTRs_to_model:n_conds_total, :]
    # print(np.max(np.abs(run_const)))
    # print(np.median((run_const)))
    
    # compute r2 of the model for each voxel
    yhat = X @ b

    ssres = np.sum((yhat - data)**2, axis=0)
    sstot = np.sum((data - np.mean(data, axis=0, keepdims=True))**2, axis=0)

    vox_r2 = 1 - (ssres/sstot)
    
    
    return hrfs, b, vox_r2
import numpy as np
import os, sys
import pandas as pd
import time
import copy

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import stats_utils
from code_utils import ordinal_regression

# another implementation that might be better
import sklearn
print(sklearn.__version__)
from code_utils.OrdinalClassifier import ordinal

def fit_withintask(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    if debug:
        subjects = [1]
    else:
        subjects = np.arange(1,8)
    n_subjects = len(subjects)
    make_time_resolved=False
    use_bigIPS = True; 
    concat_IPS = True;
    
    # first load all data for all subjects, all tasks
    maindat_all = []; repdat_all = []
    mainlabs_all = []; replabs_all = []

    for si, ss in enumerate(subjects):
       
        print('loading S%02d, main task'%ss)
        main_data, _, main_labels, roi_names = data_utils.load_main_task_data(ss, make_time_resolved, \
                                                                             use_bigIPS, concat_IPS)
        n_rois = len(roi_names)
        for ri in range(n_rois):
            # subtract mean across voxels each trial
            main_data[ri] -= np.tile(np.mean(main_data[ri], axis=1, keepdims=True), [1, main_data[ri].shape[1]])

        maindat_all += [main_data]
        mainlabs_all += [main_labels]
        
        print('loading S%02d, repeat task'%ss)
        rep_data, _, rep_labels, roi_names = data_utils.load_repeat_task_data(ss, make_time_resolved, \
                                                                             use_bigIPS, concat_IPS)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            rep_data[ri] -= np.tile(np.mean(rep_data[ri], axis=1, keepdims=True), [1, rep_data[ri].shape[1]])

        repdat_all += [rep_data]
        replabs_all += [rep_labels]

    print(roi_names)
    
    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    n_tasks = 4
    # last dimension is [x, y]
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, 2))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, 2))

    # store the predictions for individual trials
    preds_all = dict()
    probs_all = dict()
    
    for si, ss in enumerate(subjects):
        
        if debug and (si>0):
            continue
     
        preds_all[si] = dict()
        probs_all[si] = dict()
        
        # gathering labels for main task and for repeat task.
        main_labels = mainlabs_all[si]
        rep_labels = replabs_all[si]

        # all labels will be concatenated [main; repeat]
        main_grid_main = (main_labels['is_main_grid']==True) 
        main_grid_rep = (rep_labels['is_main_grid']==True) 
        is_main_grid = np.concatenate([main_grid_main, main_grid_rep], axis=0)

        inds_use_main = np.ones(np.shape(main_grid_main), dtype=bool)
        inds_use_rep= np.ones(np.shape(main_grid_rep), dtype=bool)

        xlabs_main = np.array(main_labels['ptx'])[inds_use_main]
        ylabs_main = np.array(main_labels['pty'])[inds_use_main]
        xlabs_rep = np.array(rep_labels['ptx'])[inds_use_rep]
        ylabs_rep = np.array(rep_labels['pty'])[inds_use_rep]

        xlabs = np.concatenate([xlabs_main, xlabs_rep], axis=0)
        ylabs = np.concatenate([ylabs_main, ylabs_rep], axis=0)

        # cross-validation labels, leave-one-run-out
        cv_labs_main = np.array(main_labels['run_overall'])[inds_use_main]
        cv_labs_rep = np.array(rep_labels['run_overall'])[inds_use_rep]
        cv_labs_rep += np.max(cv_labs_main)

        cv_labs = np.concatenate([cv_labs_main, cv_labs_rep], axis=0)
        n_cv = len(np.unique(cv_labs))

        # repeat task is task "4" out of 4 here
        task_labs_main = np.array(main_labels['task'])[inds_use_main]
        task_labs_rep = 4 * np.ones((np.sum(inds_use_rep), ), dtype=int)
        task_labs = np.concatenate([task_labs_main, task_labs_rep], axis=0)
        
        for ri in range(n_rois):
            
            if debug and (ri>0):
                continue
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)

            preds_all[si][ri] = dict()
            probs_all[si][ri] = dict()

            for ti, tt in enumerate([1,2,3,4]):
                
                preds_all[si][ri][ti] = dict()
                probs_all[si][ri][ti] = dict()

                tinds = task_labs==tt
                
                # loop over decoding of the x-coordinate and the y-coordinate
                # ignoring values of other coordinate
                for ii, labs in enumerate([xlabs, ylabs]):

                    labs_task = labs[tinds]
                    cv_labs_task = cv_labs[tinds]
                    is_main_grid_task = is_main_grid[tinds]

                    # data for this ROI
                    data_task = data[tinds,:]
                    dat = data_task

                    print(' processing task %d, dimension %d: %d total trials'%(tt, ii, dat.shape[0]))

                    # hold the predicted labels for entire dataset
                    nt = len(labs_task)
                    pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    prob_each = np.full(fill_value=np.nan, shape=[nt,4])

                    for cvi, cv in enumerate(np.unique(cv_labs_task)):

                        print('    fold %d of %d'%(cvi, len(np.unique(cv_labs_task))))
                        
                        # holding out one run at a time as a test set
                        # training set is all the other runs, but ONLY main grid trials.
                        trninds = (cv_labs_task!=cv) & (is_main_grid_task)
                        tstinds = cv_labs_task==cv

                        trndat = dat[trninds,:]
                        tstdat = dat[tstinds,:]

                        trnlabs = (labs_task[trninds]*10).astype(int)
                        assert(len(np.unique(trnlabs))==4)

                        C = 0.1;
                        model_pars = [C, n_threads]
                        m = ordinal_regression.ordinal_regress_model(ordinal_regression.get_model)
                        m.fit(trndat, trnlabs, model_pars)
                        
                        # C = 0.1;
                        # model_pars = [C, n_threads]
                        # model = ordinal_regression.get_model(*model_pars)
                        # m = ordinal.OrdinalClassifier(model)
                        # m.fit(trndat, trnlabs)

                        pred, prob = m.predict(tstdat)
                        pred = pred/10
                        
                        print(prob[0,:])
                        print(prob.shape)
                        assert(np.all(np.sum(prob, axis=1).round(9)==1))
                        
                        # pred is the categorical prediction, prob is continuous
                        pred_labs[tstinds] = pred 
                        prob_each[tstinds] = prob
                        
                    if not debug:
                        assert(not np.any(np.isnan(pred_labs)))
                        assert(not np.any(np.isnan(prob_each)))
                        
                    acc_bytask[si,ri,ti,ii] = np.mean(pred_labs[is_main_grid_task]==\
                                                                        labs_task[is_main_grid_task])
                    dprime_bytask[si,ri,ti,ii] = stats_utils.get_dprime(pred_labs[is_main_grid_task], \
                                                                        labs_task[is_main_grid_task])
                    
                    print('acc is %.2f'%(acc_bytask[si,ri,ti,ii]))
                    
                    # save trial-wise predictions and probability scores
                    preds_all[si][ri][ti][ii] = pred_labs
                    probs_all[si][ri][ti][ii] = prob_each
                    
                    
        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        save_filename = os.path.join(save_folder, 'ordinal_regression_withintask.npy')
        # save_filename = os.path.join(save_folder, 'ordinal_regression_method2_withintask.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'roi_names': roi_names, \
                               })
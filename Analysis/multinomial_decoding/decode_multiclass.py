import numpy as np
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis, sklearn.linear_model
from joblib import effective_n_jobs
import time
import copy
import datetime

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import decoding_utils
from code_utils import stats_utils

def decode_withintask(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    print('cpus available = %d'%(effective_n_jobs(-1)))
    
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
    # penalties to eval
    c_values = np.logspace(-9, 1, 20)
    n_grid_pts = 16
    
    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    n_tasks = 4
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks))

    n_cv = 12;
    # store which c value is best
    # these accs are on the nested training data, not the testing data
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks, n_cv), np.nan)                        
   
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
        
        pt_labs = np.array([xlabs, ylabs]).T
        grid_pts, grid_labs_main, counts = np.unique(pt_labs[is_main_grid], axis=0, return_inverse=True, return_counts=True)
        assert(n_grid_pts==grid_pts.shape[0])
        assert(np.all(counts==counts[0]))
        # grid labs are nan for the off-grid trials
        grid_labs_all = np.full(fill_value=np.nan, shape=np.shape(is_main_grid))
        grid_labs_all[is_main_grid] = grid_labs_main

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
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>3):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            preds_all[si][ri] = dict()
            probs_all[si][ri] = dict()

            for ti, tt in enumerate([1,2,3,4]):
                
                tinds = task_labs==tt
                
                grid_labs_task = grid_labs_all[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                print(' processing task %d: %d total trials'%(tt, dat.shape[0]))
                
                # hold the predicted labels for this task
                nt = len(grid_labs_task)
                pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                prob_each = np.full(fill_value=np.nan, shape=[nt,n_grid_pts])

                for cvi, cv in enumerate(np.unique(cv_labs_task)):

                    if debug & (cvi>0):
                        continue
                    
                    # holding out one run at a time as a test set
                    # training set is all the other runs, only main grid trials.
                    trninds = (cv_labs_task!=cv) & (is_main_grid_task)
                    tstinds = cv_labs_task==cv

                    trndat = dat[trninds,:]
                    tstdat = dat[tstinds,:]
                    
                    trnlabs = grid_labs_task[trninds]
                    assert(not np.any(np.isnan(trnlabs)))
                    
                    print(trndat.shape, tstdat.shape)
                    
                    # # do regularization parameter (c) selection
                    # # this is based on training data only, for the current fold.
                    # # cross-validate using leave-one-run-out (for just the training runs here)
                    nest_cv_labs = cv_labs_task[trninds]
                    # n_cv_nest = len(np.unique(nest_cv_labs))
                    nest_cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                    nest_cv_generator = nest_cv_obj.split(trndat, trnlabs, nest_cv_labs)
                    # n_cv_nest = 10
                    
                    # define model
                    st = time.time()
                    model = sklearn.linear_model.LogisticRegressionCV(\
                                                                      # cv = n_cv_nest, \
                                                                      cv = nest_cv_generator, \
                                                                    Cs = c_values, \
                                                                    multi_class='multinomial',\
                                                                    solver='lbfgs', \
                                                                    penalty='l2', \
                                                                    n_jobs = n_threads , \
                                                                    max_iter = 1000)
                    model.fit(trndat, trnlabs)
                    elapsed = time.time() - st

                    # pull out the accuracy of the model for each C value
                    # averaging across the nested CV folds
                    a = np.mean(model.scores_[0], axis=0)
                    c = model.C_[0]
                    assert(c_values[np.argmax(a)]==c)

                    print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                    sys.stdout.flush()

                    acc_each_cval[si, ri, ti, cvi,:] = a
                    best_cval[si, ri, ti, cvi] = c

                    # finally, predict on the held-out test data here
                    pred = model.predict(tstdat)
                    prob = model.predict_proba(tstdat)
                    # pred is the categorical prediction, prob is continuous
                    assert(np.all(np.sum(prob, axis=1).round(9)==1))
                    assert(np.all(np.argmax(prob,axis=1)==pred.astype(int)))

                    pred_labs[tstinds] = pred
                    prob_each[tstinds,:] = prob

                if not debug:
                    assert(not np.any(np.isnan(pred_labs)))
                    assert(not np.any(np.isnan(prob_each)))
                
                # save trial-wise predictions and probability scores
                preds_all[si][ri][ti] = pred_labs
                probs_all[si][ri][ti] = prob_each
                
                # compute some performance metrics
                # these are just the "main grid" because easy to compute accuracy
                assert(not np.any(np.isnan(grid_labs_task[is_main_grid_task])))
                acc_bytask[si,ri,ti] = np.mean(pred_labs[is_main_grid_task]==grid_labs_task[is_main_grid_task])
                dprime_bytask[si,ri,ti] = stats_utils.get_dprime(pred_labs[is_main_grid_task], grid_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_multiclass_withintask_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_multiclass_withintask.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'grid_pts': grid_pts, \
                               'roi_names': roi_names, \
                               })


        
        
def decode_withintask_permutationtest(debug=False, n_threads=8, n_iter=1000, rndseed=None):
    
    if rndseed is None:
        rndseed = int(datetime.datetime.now().strftime('%f'))
        
    print('debug = %s, n_threads = %d, n_iter = %d, rndseed = %d'%(debug, n_threads, n_iter, rndseed))
    print('cpus available = %d'%(effective_n_jobs(-1)))
    
    np.random.seed(rndseed)
    
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
    # penalties to eval
    # c_values = np.logspace(-9, 1, 20)
    c_value = 0.023 
    # this value is the median of values that were best for real data
    # it's too slow to pick a value for every shuffle iteration
    
    n_grid_pts = 16
    
    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    n_tasks = 4
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_iter))
    # dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_iter))

    st_allsubs = time.time()
    
    for si, ss in enumerate(subjects):
        
        if debug and (si>0):
            continue
     
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
        
        pt_labs = np.array([xlabs, ylabs]).T
        grid_pts, grid_labs_main, counts = np.unique(pt_labs[is_main_grid], axis=0, return_inverse=True, return_counts=True)
        assert(n_grid_pts==grid_pts.shape[0])
        assert(np.all(counts==counts[0]))
        # grid labs are nan for the off-grid trials
        grid_labs_all = np.full(fill_value=np.nan, shape=np.shape(is_main_grid))
        grid_labs_all[is_main_grid] = grid_labs_main

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
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>3):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            
            for ti, tt in enumerate([1,2,3,4]):
                
                st_overall = time.time()
                
                tinds = (task_labs==tt) & is_main_grid
                
                grid_labs_task = grid_labs_all[tinds]
                cv_labs_task = cv_labs[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                
                print('Processing task %d: %d total trials'%(tt, dat.shape[0]))
                
                # hold the predicted labels for this task
                nt = len(grid_labs_task)
                
                # make shuffled labels for each iteration of perm test
                # shuffle within task runs
                shuff_labs_all = np.full(shape=(nt,n_iter), fill_value=np.nan)
                for cvi, cv in enumerate(np.unique(cv_labs_task)):
                    
                    # shuffle the "main grid" labels only here
                    # classifier training set is only main grid
                    # the not-main grid trials don't get used at all in this code
                    # cvidx1 = (cv_labs_task==cv) & is_main_grid_task
                    cvidx1 = (cv_labs_task==cv)
                    assert(np.all(np.isnan(shuff_labs_all[cvidx1])))
                    actual_labs1 = grid_labs_task[cvidx1]
                    
                    for xi in range(n_iter):
                        shuff_labs_all[cvidx1,xi] = actual_labs1[np.random.permutation(len(actual_labs1))]
                        
                # assert(not np.any(np.isnan(shuff_labs_all[is_main_grid_task,:])))
                assert(not np.any(np.isnan(shuff_labs_all)))
                
                
                for xi in range(n_iter):
                    
                    print('shuffle iter %d of %d'%(xi, n_iter))
                          
                    shuff_labs = shuff_labs_all[:,xi]
                    
                    pred_labs_shuff = np.full(fill_value=np.nan, shape=[nt,])
                
                    # to cross-validate, using the built in cross-validator
                    # from sklearn, and pass this into the LogisticRegressionCV funct.
                    # this is faster than doing it manually
                    cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                    cv_generator = cv_obj.split(dat, grid_labs_task, cv_labs_task)
                        
                    # define model
                    st = time.time()
                    model = sklearn.linear_model.LogisticRegressionCV(cv = cv_generator, \
                                                                      Cs = [c_value], \
                                                                    multi_class='multinomial',\
                                                                    solver='lbfgs', \
                                                                    penalty='l2', \
                                                                    n_jobs = n_threads , \
                                                                    max_iter = 1000, \
                                                                    refit = False)
                    model.fit(dat, shuff_labs)
                    elapsed = time.time() - st
                    
                    # the accuracy for each fold is in .scores_
                    # this is cross-validated, so can use it directly 
                    # since there is only one C value this is ok
                    acc_each_cv = model.scores_[0]
                    acc = np.mean(acc_each_cv[:,0])
                   
                    acc_bytask[si,ri,ti,xi] = acc
        
                et_overall = time.time()
                elapsed = et_overall - st_overall
                print('Took %.6f s for S%02d %s task %d (average acc = %.2f)'%(elapsed, ss, roi_names[ri], \
                                                                               tt, np.mean(acc_bytask[si,ri,ti,:])))
            
        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_multiclass_withintask_permutationtest_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_multiclass_withintask_permutationtest.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               # 'dprime_bytask': dprime_bytask, \
                               'grid_pts': grid_pts, \
                               'roi_names': roi_names, \
                               'rndseed': rndseed, \
                               })

    et_allsubs = time.time()
    elapsed = et_allsubs - st_allsubs
    print('Took %.6f s for all subs'%(elapsed))
        
       
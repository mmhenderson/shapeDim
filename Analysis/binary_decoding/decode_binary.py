import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis, sklearn.linear_model
import time
from joblib import effective_n_jobs
import datetime
import scipy.stats


# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import decoding_utils
from code_utils import stats_utils


def decode_withintask(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    
    subjects = np.arange(1,11)
        
    if debug:
        subjects = subjects[0:1]
        
    print(subjects)
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
    
    n_tasks = 4
    
    # three different ways to do binary decoding
    n_bounds = 3;
    bound_names = ['Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'];
    quad_groups = [[[1, 4], [2, 3]],
                    [[1, 2], [3, 4]],
                    [[1, 3], [2, 4]]];

    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))

    n_cv = 12;
    # store which c value is best
    # these accs are on the nested training data, not the testing data
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv), np.nan)                        
   
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
        concat_labels = pd.concat([main_labels, rep_labels], axis=0)
        
        # cross-validation labels, leave-one-run-out
        cv_labs_main = np.array(main_labels['run_overall'])
        cv_labs_rep = np.array(rep_labels['run_overall'])
        cv_labs_rep += np.max(cv_labs_main)
        
        cv_labs = np.concatenate([cv_labs_main, cv_labs_rep], axis=0)
        
        is_main_grid = np.array(concat_labels['is_main_grid']==True) 
        quad_labs = np.array(concat_labels['quadrant'])
        
        task_labs = np.array(concat_labels['task'])
        # repeat task is task "4" out of 4 here
        task_labs[~np.isin(task_labs, [1,2,3])] = 4
        
        
        for ri in range(n_rois):
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            preds_all[si][ri] = dict()
            probs_all[si][ri] = dict()

            # loop over tasks
            for ti, tt in enumerate([1,2,3,4]):
                
                preds_all[si][ri][ti] = dict()
                probs_all[si][ri][ti] = dict()
                
                tinds = task_labs==tt
                
                quad_labs_task = quad_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                # loop over ways of doing decoding, within each task
                for bi in range(n_bounds):

                    # convert quadrants to binary labels, for the current boundary
                    categ_labs_task = np.zeros(np.shape(quad_labs_task))
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][0])] = 1
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][1])] = 2
                    assert(not np.any(categ_labs_task==0))

                    print(' processing task %d, boundary %d: %d total trials'%(tt, bi, dat.shape[0]))

                    # hold the predicted labels for this task
                    nt = len(categ_labs_task)
                    pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    prob_each = np.full(fill_value=np.nan, shape=[nt,2])

                    for cvi, cv in enumerate(np.unique(cv_labs_task)):

                        if debug & (cvi>0):
                            continue

                        # holding out one run at a time as a test set
                        # training set is all the other runs, only main grid trials.
                        trninds = (cv_labs_task!=cv) & (is_main_grid_task)
                        tstinds = cv_labs_task==cv

                        trndat = dat[trninds,:]
                        tstdat = dat[tstinds,:]

                        trnlabs = categ_labs_task[trninds]
                        assert(not np.any(np.isnan(trnlabs)))
                        # check for balanced labels
                        un, counts = np.unique(categ_labs_task[trninds], return_counts=True)
                        assert(counts[0]==counts[1])

                        print(trndat.shape, tstdat.shape)

                        # do regularization parameter (c) selection
                        # this is based on training data only, for the current fold.
                        # cross-validate using leave-one-run-out (for just the training runs here)
                        nest_cv_labs = cv_labs_task[trninds]
                        nest_cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                        nest_cv_generator = nest_cv_obj.split(trndat, trnlabs, nest_cv_labs)

                        # define model
                        st = time.time()
                        model = sklearn.linear_model.LogisticRegressionCV(cv = nest_cv_generator, \
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
                        # print(np.unique(trnlabs))
                        # print(model.scores_.keys())
                        # print(model.scores_)
                        # print(model.scores_.shape)
                        a = np.mean(model.scores_[un[1]], axis=0)
                        c = model.C_
                        assert(c_values[np.argmax(a)]==c)

                        print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                        sys.stdout.flush()

                        acc_each_cval[si, ri, ti, bi, cvi,:] = a
                        best_cval[si, ri, ti, bi, cvi] = c

                        # finally, predict on the held-out test data here
                        pred = model.predict(tstdat)
                        prob = model.predict_proba(tstdat)
                        # pred is the categorical prediction, prob is continuous
                        assert(np.all(np.sum(prob, axis=1).round(9)==1))
                        assert(np.all(un[np.argmax(prob,axis=1)]==pred))
                        
                        pred_labs[tstinds] = pred
                        prob_each[tstinds,:] = prob

                    if not debug:
                        assert(not np.any(np.isnan(pred_labs)))
                        assert(not np.any(np.isnan(prob_each)))

                    # save trial-wise predictions and probability scores
                    preds_all[si][ri][ti][bi] = pred_labs
                    probs_all[si][ri][ti][bi] = prob_each

                    # compute some performance metrics
                    # these are just the "main grid" because easy to compute accuracy
                    assert(not np.any(np.isnan(categ_labs_task[is_main_grid_task])))
                    acc_bytask[si,ri,ti,bi] = \
                        np.mean(pred_labs[is_main_grid_task]==categ_labs_task[is_main_grid_task])
                    dprime_bytask[si,ri,ti,bi] = \
                        stats_utils.get_dprime(pred_labs[is_main_grid_task], categ_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'roi_names': roi_names, \
                               })


        
def decode_withintask_permutationtest(debug=False, n_threads=8, n_iter=1000, rndseed=None,):
    
    if rndseed is None:
        rndseed = int(datetime.datetime.now().strftime('%f'))
        
    print('debug = %s, n_threads = %d, n_iter = %d, rndseed = %d'%(debug, n_threads, n_iter, rndseed))
    print('cpus available = %d'%(effective_n_jobs(-1)))
    
    np.random.seed(rndseed)
    
    subjects = np.arange(1,11)
        
    if debug:
        subjects = subjects[0:1]
        
    print(subjects)
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
    c_value = 0.007
    # this value is the median of values that were best for real data
    # it's too slow to pick a value for every shuffle iteration
    
    n_tasks = 4
    
    # three different ways to do binary decoding
    n_bounds = 3;
    bound_names = ['Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'];
    quad_groups = [[[1, 4], [2, 3]],
                    [[1, 2], [3, 4]],
                    [[1, 3], [2, 4]]];

    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds, n_iter))
    # dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))

    # n_cv = 12;
    # # store which c value is best
    # # these accs are on the nested training data, not the testing data
    # acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv, len(c_values)), np.nan)
    # best_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv), np.nan)                        
   
    # # store the predictions for individual trials
    # preds_all = dict()
    # probs_all = dict()

    for si, ss in enumerate(subjects):
        
        if debug and (si>0):
            continue
            
#         preds_all[si] = dict()
#         probs_all[si] = dict()
        
        # gathering labels for main task and for repeat task.
        main_labels = mainlabs_all[si]
        rep_labels = replabs_all[si]
        # all labels will be concatenated [main; repeat]
        concat_labels = pd.concat([main_labels, rep_labels], axis=0)
        
        # cross-validation labels, leave-one-run-out
        cv_labs_main = np.array(main_labels['run_overall'])
        cv_labs_rep = np.array(rep_labels['run_overall'])
        cv_labs_rep += np.max(cv_labs_main)
        
        cv_labs = np.concatenate([cv_labs_main, cv_labs_rep], axis=0)
        
        is_main_grid = np.array(concat_labels['is_main_grid']==True) 
        quad_labs = np.array(concat_labels['quadrant'])
        
        task_labs = np.array(concat_labels['task'])
        # repeat task is task "4" out of 4 here
        task_labs[~np.isin(task_labs, [1,2,3])] = 4
        
        
        for ri in range(n_rois):
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
#             preds_all[si][ri] = dict()
#             probs_all[si][ri] = dict()

            # loop over tasks
            for ti, tt in enumerate([1,2,3,4]):
                
#                 preds_all[si][ri][ti] = dict()
#                 probs_all[si][ri][ti] = dict()
                
                tinds = task_labs==tt
                
                quad_labs_task = quad_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                # loop over ways of doing decoding, within each task
                for bi in range(n_bounds):
                    
                    st_overall = time.time()
                

                    # convert quadrants to binary labels, for the current boundary
                    categ_labs_task = np.zeros(np.shape(quad_labs_task))
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][0])] = 1
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][1])] = 2
                    assert(not np.any(categ_labs_task==0))

                    print(' processing task %d, boundary %d: %d total trials'%(tt, bi, dat.shape[0]))

                    # hold the predicted labels for this task
                    nt = len(categ_labs_task)
                    # pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    # prob_each = np.full(fill_value=np.nan, shape=[nt,2])

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
                        actual_labs1 = categ_labs_task[cvidx1]

                        for xi in range(n_iter):
                            shuff_labs_all[cvidx1,xi] = actual_labs1[np.random.permutation(len(actual_labs1))]
                    
                    assert(not np.any(np.isnan(shuff_labs_all)))
                
                    for xi in range(n_iter):
                    
                        print('shuffle iter %d of %d'%(xi, n_iter))

                        shuff_labs = shuff_labs_all[:,xi]

                        # pred_labs_shuff = np.full(fill_value=np.nan, shape=[nt,])

                        # to cross-validate, using the built in cross-validator
                        # from sklearn, and pass this into the LogisticRegressionCV funct.
                        # this is faster than doing it manually
                        cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                        cv_generator = cv_obj.split(dat, categ_labs_task, cv_labs_task)

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
                        
                        # print(np.unique(shuff_labs))
                        # the accuracy for each fold is in .scores_
                        # this is cross-validated, so can use it directly 
                        # since there is only one C value this is ok
                        keys = list(model.scores_.keys())
                        assert(len(keys)==1)
                        key = keys[0]
                        # print(model.scores_)
                        # print(model.scores_.keys())
                        # print(model.scores_[key].shape)
                        
                        acc_each_cv = model.scores_[key]
                        acc = np.mean(acc_each_cv[:,0])

                        acc_bytask[si,ri,ti,bi,xi] = acc
                        
                    et_overall = time.time()
                    elapsed = et_overall - st_overall
                    print('Took %.6f s for S%02d %s task %d bound %d (average acc = %.2f)'%(elapsed, ss, roi_names[ri], \
                                                                               tt, bi, np.mean(acc_bytask[si,ri,ti,bi,:])))
            

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_permutationtest_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_permutationtest.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'roi_names': roi_names, \
                               'rndseed': rndseed, \
                               })

        
        

def decode_withintask_match_num_voxels(debug=False, n_threads=8, rndseed=None):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    if rndseed is None:
        rndseed = int(datetime.datetime.now().strftime('%f'))
        
    print('rndseed = %d'%(rndseed))
    
    np.random.seed(rndseed)
    
    subjects = np.arange(1,11)
        
    if debug:
        subjects = subjects[0:1]
        
    print(subjects)
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
        
        
    # I'm going to use same num voxels across all ROIs. 
    # Doing this matching within a subject, so that we can still maximize power.
    # Use data from the silhouette localizer task to get t-statistics.

    vox_use_all = []

    for si, ss in enumerate(subjects):

        n_vox_each_roi = np.array([m.shape[1] for m in repdat_all[si]])
        min_vox_this_sub = np.min(n_vox_each_roi)

        print('\nFor S%02d, the minimum number of voxels: is %d'%(ss, min_vox_this_sub))

        v = []

        tvals = get_tstats(ss)

        for ri in range(n_rois):

            print('%s: choosing the best %d voxels to use'%(roi_names[ri], min_vox_this_sub))

            t = tvals[ri][:,0]

            best_vox_order = np.flip(np.argsort(t))
            vox_use = best_vox_order[0:min_vox_this_sub]

            v += [vox_use]

            # printing max and min. make sure this selection worked
            print(np.min(t[vox_use]), np.max(t[vox_use]))

        vox_use_all += [v]


    print(roi_names)
    # penalties to eval
    c_values = np.logspace(-9, 1, 20)
    
    n_tasks = 4
    
    # three different ways to do binary decoding
    n_bounds = 3;
    bound_names = ['Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'];
    quad_groups = [[[1, 4], [2, 3]],
                    [[1, 2], [3, 4]],
                    [[1, 3], [2, 4]]];

    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))

    n_cv = 12;
    # store which c value is best
    # these accs are on the nested training data, not the testing data
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv), np.nan)                        
   
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
        concat_labels = pd.concat([main_labels, rep_labels], axis=0)
        
        # cross-validation labels, leave-one-run-out
        cv_labs_main = np.array(main_labels['run_overall'])
        cv_labs_rep = np.array(rep_labels['run_overall'])
        cv_labs_rep += np.max(cv_labs_main)
        
        cv_labs = np.concatenate([cv_labs_main, cv_labs_rep], axis=0)
        
        is_main_grid = np.array(concat_labels['is_main_grid']==True) 
        quad_labs = np.array(concat_labels['quadrant'])
        
        task_labs = np.array(concat_labels['task'])
        # repeat task is task "4" out of 4 here
        task_labs[~np.isin(task_labs, [1,2,3])] = 4
        
        
        for ri in range(n_rois):
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            preds_all[si][ri] = dict()
            probs_all[si][ri] = dict()

            # loop over tasks
            for ti, tt in enumerate([1,2,3,4]):
                
                preds_all[si][ri][ti] = dict()
                probs_all[si][ri][ti] = dict()
                
                tinds = task_labs==tt
                
                quad_labs_task = quad_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                # loop over ways of doing decoding, within each task
                for bi in range(n_bounds):

                    # convert quadrants to binary labels, for the current boundary
                    categ_labs_task = np.zeros(np.shape(quad_labs_task))
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][0])] = 1
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][1])] = 2
                    assert(not np.any(categ_labs_task==0))

                    print(' processing task %d, boundary %d: %d total trials'%(tt, bi, dat.shape[0]))

                    # hold the predicted labels for this task
                    nt = len(categ_labs_task)
                    pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    prob_each = np.full(fill_value=np.nan, shape=[nt,2])

                    for cvi, cv in enumerate(np.unique(cv_labs_task)):

                        if debug & (cvi>0):
                            continue

                        # holding out one run at a time as a test set
                        # training set is all the other runs, only main grid trials.
                        trninds = (cv_labs_task!=cv) & (is_main_grid_task)
                        tstinds = cv_labs_task==cv

                        trndat = dat[trninds,:]
                        tstdat = dat[tstinds,:]

                        trnlabs = categ_labs_task[trninds]
                        assert(not np.any(np.isnan(trnlabs)))
                        # check for balanced labels
                        un, counts = np.unique(categ_labs_task[trninds], return_counts=True)
                        assert(counts[0]==counts[1])

                        print(trndat.shape, tstdat.shape)
                        
                        vox_use = vox_use_all[si][ri]
                        
                        print('choosing the best %d voxels to use'%len(vox_use))
                        trndat = trndat[:,vox_use]
                        tstdat = tstdat[:,vox_use]
                        
                        print(trndat.shape, tstdat.shape)


                        # do regularization parameter (c) selection
                        # this is based on training data only, for the current fold.
                        # cross-validate using leave-one-run-out (for just the training runs here)
                        nest_cv_labs = cv_labs_task[trninds]
                        nest_cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                        nest_cv_generator = nest_cv_obj.split(trndat, trnlabs, nest_cv_labs)

                        # define model
                        st = time.time()
                        model = sklearn.linear_model.LogisticRegressionCV(cv = nest_cv_generator, \
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
                        # print(np.unique(trnlabs))
                        # print(model.scores_.keys())
                        # print(model.scores_)
                        # print(model.scores_.shape)
                        a = np.mean(model.scores_[un[1]], axis=0)
                        c = model.C_
                        assert(c_values[np.argmax(a)]==c)

                        print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                        sys.stdout.flush()

                        acc_each_cval[si, ri, ti, bi, cvi,:] = a
                        best_cval[si, ri, ti, bi, cvi] = c

                        # finally, predict on the held-out test data here
                        pred = model.predict(tstdat)
                        prob = model.predict_proba(tstdat)
                        # pred is the categorical prediction, prob is continuous
                        assert(np.all(np.sum(prob, axis=1).round(9)==1))
                        assert(np.all(un[np.argmax(prob,axis=1)]==pred))
                        
                        pred_labs[tstinds] = pred
                        prob_each[tstinds,:] = prob

                    if not debug:
                        assert(not np.any(np.isnan(pred_labs)))
                        assert(not np.any(np.isnan(prob_each)))

                    # save trial-wise predictions and probability scores
                    preds_all[si][ri][ti][bi] = pred_labs
                    probs_all[si][ri][ti][bi] = prob_each

                    # compute some performance metrics
                    # these are just the "main grid" because easy to compute accuracy
                    assert(not np.any(np.isnan(categ_labs_task[is_main_grid_task])))
                    acc_bytask[si,ri,ti,bi] = \
                        np.mean(pred_labs[is_main_grid_task]==categ_labs_task[is_main_grid_task])
                    dprime_bytask[si,ri,ti,bi] = \
                        stats_utils.get_dprime(pred_labs[is_main_grid_task], categ_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_matchnumvoxels_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_matchnumvoxels.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'roi_names': roi_names, \
                               'rndseed': rndseed,
                               })


def get_tstats(ss):
    
    # this is from the sihouette localizer task - these are the t-statistics
    # for a GLM in that task. Same data that is used to create our ROIs by thresholding.
    # Note this isn't used to threshold IPS, so IPS is bigger and has negative t-values in it.
    
    # i am only using this for the "match num voxels" functions
    
    # these came from GetLocVals_bigIPS.m
    locvals_fn = os.path.join(root, 'Samples','LocVals_bigIPS_S%02d.mat'%ss)
    print(locvals_fn)
    loc = file_utils.load_locvalsfile_h5py(locvals_fn)
    
    roi_names = loc['ROI_names']
    n_rois = len(roi_names)
    n_hemis = len(loc['hemis'])
    
    tvals = []

    for rr in range(n_rois):
        
        print(roi_names[rr])

        tvals_this_roi = []

        for hh in range(n_hemis):

            roi_inds_num = np.array(loc['ROIs']['voxel_inds'][rr][hh])

            inds_this_hemi = np.isin(np.array(loc['all_vox_concat']), roi_inds_num)[:,0]

            if np.sum(inds_this_hemi)>0:

                tvals_this_hemi = loc['loc_vals'][inds_this_hemi]
                tvals_this_roi += [tvals_this_hemi]

        tvals_this_roi = np.concatenate(tvals_this_roi, axis=0)
        print(tvals_this_roi.shape)

        print(np.min(tvals_this_roi), np.max(tvals_this_roi))
        nVox = tvals_this_roi.shape[1]
        
        tvals += [tvals_this_roi]
        
    # now concatenating the IPS subregions
    # this is what is also done in data_utils.load_main_task_data
    tvals_adj = tvals[0:5] + tvals[9:]
    ips_roi = [5,6,7,8]
    ips_concat = np.concatenate([tvals[rr] for rr in ips_roi], axis=0)
    ips_concat.shape
    
    tvals_adj += [ips_concat]
    
    print('final sizes:')
    print([t.shape[0] for t in tvals_adj])
        
    return tvals_adj


        

def decode_withintask_match_num_voxels_method2(debug=False, n_threads=8, rndseed=None):
    
    # this is a different method, using repeat task to choose voxels
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    if rndseed is None:
        rndseed = int(datetime.datetime.now().strftime('%f'))
        
    print('rndseed = %d'%(rndseed))
    
    np.random.seed(rndseed)
    
    subjects = np.arange(1,11)
        
    if debug:
        subjects = subjects[0:1]
        
    print(subjects)
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
        
        
    # I'm going to use same num voxels across all ROIs. 
    # Doing this matching within a subject, so that we can still maximize power.
    # Use data from the "repeat" task to choose voxels, so it's independent 
    # of the categorization tasks.
   
    vox_use_all = []

    for si, ss in enumerate(subjects):
        
        n_vox_each_roi = np.array([m.shape[1] for m in repdat_all[si]])
        min_vox_this_sub = np.min(n_vox_each_roi)

        print('\nFor S%02d, the minimum number of voxels: is %d'%(ss, min_vox_this_sub))

        v = []

        for ri in range(n_rois):

            print('%s: choosing the best %d voxels to use'%(roi_names[ri], min_vox_this_sub))

            # Controlling the number of voxels in each ROI, per subject
            # To choose voxels, I'm doing a one-way anova on the repeat task data.
            # Use "quadrant" as the labels for this.
            inds_use = np.array(replabs_all[si]['is_main_grid'])==1
            q = np.array(replabs_all[si]['quadrant'])[inds_use]

            d = repdat_all[si][ri][inds_use,:]

            samples_for_anova = [d[q==ii,:] for ii in [1,2,3,4]]

            f_result = scipy.stats.f_oneway(*samples_for_anova, axis=0)
            voxel_f = f_result.statistic

            best_vox_order = np.flip(np.argsort(voxel_f))
            vox_use = best_vox_order[0:min_vox_this_sub]

            v += [vox_use]

        vox_use_all += [v]


    print(roi_names)
    # penalties to eval
    c_values = np.logspace(-9, 1, 20)
    
    n_tasks = 4
    
    # three different ways to do binary decoding
    n_bounds = 3;
    bound_names = ['Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'];
    quad_groups = [[[1, 4], [2, 3]],
                    [[1, 2], [3, 4]],
                    [[1, 3], [2, 4]]];

    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_bounds))

    n_cv = 12;
    # store which c value is best
    # these accs are on the nested training data, not the testing data
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks, n_bounds, n_cv), np.nan)                        
   
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
        concat_labels = pd.concat([main_labels, rep_labels], axis=0)
        
        # cross-validation labels, leave-one-run-out
        cv_labs_main = np.array(main_labels['run_overall'])
        cv_labs_rep = np.array(rep_labels['run_overall'])
        cv_labs_rep += np.max(cv_labs_main)
        
        cv_labs = np.concatenate([cv_labs_main, cv_labs_rep], axis=0)
        
        is_main_grid = np.array(concat_labels['is_main_grid']==True) 
        quad_labs = np.array(concat_labels['quadrant'])
        
        task_labs = np.array(concat_labels['task'])
        # repeat task is task "4" out of 4 here
        task_labs[~np.isin(task_labs, [1,2,3])] = 4
        
        
        for ri in range(n_rois):
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            preds_all[si][ri] = dict()
            probs_all[si][ri] = dict()

            # loop over tasks
            for ti, tt in enumerate([1,2,3,4]):
                
                preds_all[si][ri][ti] = dict()
                probs_all[si][ri][ti] = dict()
                
                tinds = task_labs==tt
                
                quad_labs_task = quad_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                # loop over ways of doing decoding, within each task
                for bi in range(n_bounds):

                    # convert quadrants to binary labels, for the current boundary
                    categ_labs_task = np.zeros(np.shape(quad_labs_task))
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][0])] = 1
                    categ_labs_task[np.isin(quad_labs_task, quad_groups[bi][1])] = 2
                    assert(not np.any(categ_labs_task==0))

                    print(' processing task %d, boundary %d: %d total trials'%(tt, bi, dat.shape[0]))

                    # hold the predicted labels for this task
                    nt = len(categ_labs_task)
                    pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    prob_each = np.full(fill_value=np.nan, shape=[nt,2])

                    for cvi, cv in enumerate(np.unique(cv_labs_task)):

                        if debug & (cvi>0):
                            continue

                        # holding out one run at a time as a test set
                        # training set is all the other runs, only main grid trials.
                        trninds = (cv_labs_task!=cv) & (is_main_grid_task)
                        tstinds = cv_labs_task==cv

                        trndat = dat[trninds,:]
                        tstdat = dat[tstinds,:]

                        trnlabs = categ_labs_task[trninds]
                        assert(not np.any(np.isnan(trnlabs)))
                        # check for balanced labels
                        un, counts = np.unique(categ_labs_task[trninds], return_counts=True)
                        assert(counts[0]==counts[1])

                        print(trndat.shape, tstdat.shape)
                        
                        vox_use = vox_use_all[si][ri]
                        
                        print('choosing the best %d voxels to use'%len(vox_use))
                        trndat = trndat[:,vox_use]
                        tstdat = tstdat[:,vox_use]
                        
                        print(trndat.shape, tstdat.shape)


                        # do regularization parameter (c) selection
                        # this is based on training data only, for the current fold.
                        # cross-validate using leave-one-run-out (for just the training runs here)
                        nest_cv_labs = cv_labs_task[trninds]
                        nest_cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                        nest_cv_generator = nest_cv_obj.split(trndat, trnlabs, nest_cv_labs)

                        # define model
                        st = time.time()
                        model = sklearn.linear_model.LogisticRegressionCV(cv = nest_cv_generator, \
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
                        # print(np.unique(trnlabs))
                        # print(model.scores_.keys())
                        # print(model.scores_)
                        # print(model.scores_.shape)
                        a = np.mean(model.scores_[un[1]], axis=0)
                        c = model.C_
                        assert(c_values[np.argmax(a)]==c)

                        print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                        sys.stdout.flush()

                        acc_each_cval[si, ri, ti, bi, cvi,:] = a
                        best_cval[si, ri, ti, bi, cvi] = c

                        # finally, predict on the held-out test data here
                        pred = model.predict(tstdat)
                        prob = model.predict_proba(tstdat)
                        # pred is the categorical prediction, prob is continuous
                        assert(np.all(np.sum(prob, axis=1).round(9)==1))
                        assert(np.all(un[np.argmax(prob,axis=1)]==pred))
                        
                        pred_labs[tstinds] = pred
                        prob_each[tstinds,:] = prob

                    if not debug:
                        assert(not np.any(np.isnan(pred_labs)))
                        assert(not np.any(np.isnan(prob_each)))

                    # save trial-wise predictions and probability scores
                    preds_all[si][ri][ti][bi] = pred_labs
                    probs_all[si][ri][ti][bi] = prob_each

                    # compute some performance metrics
                    # these are just the "main grid" because easy to compute accuracy
                    assert(not np.any(np.isnan(categ_labs_task[is_main_grid_task])))
                    acc_bytask[si,ri,ti,bi] = \
                        np.mean(pred_labs[is_main_grid_task]==categ_labs_task[is_main_grid_task])
                    dprime_bytask[si,ri,ti,bi] = \
                        stats_utils.get_dprime(pred_labs[is_main_grid_task], categ_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if debug:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_matchnumvoxels2_DEBUG.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_binary_withintask_matchnumvoxels2.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'roi_names': roi_names, \
                               'rndseed': rndseed,
                               })


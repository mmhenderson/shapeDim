import numpy as np
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis, sklearn.linear_model
import time
import copy

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import decoding_utils
from code_utils import stats_utils

def decode_withintask(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
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
   
    pred_labs_all = dict()

    for si, ss in enumerate(subjects):
        
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

        n_total_trials = len(is_main_grid)
        
        pred_labs_all[ss] = np.zeros((n_rois, n_total_trials))
        
        for ri in range(n_rois):
            
            main_data = maindat_all[si][ri]
            rep_data = repdat_all[si][ri]
            data = np.concatenate([main_data, rep_data], axis=0)
            print(data.shape)
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            for ti, tt in enumerate([1,2,3,4]):
                
                tinds = task_labs==tt
                
                grid_labs_task = grid_labs_all[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                print(' processing task %d: %d total trials'%(tt, dat.shape[0]))
                
                # hold the predicted labels for entire dataset
                pred_labs = np.full(fill_value=np.nan, shape=(np.shape(grid_labs_task)))

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
                    a = np.mean(model.scores_[0], axis=0)
                    c = model.C_[0]
                    assert(c_values[np.argmax(a)]==c)

                    print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                    sys.stdout.flush()

                    acc_each_cval[si, ri, ti, cvi,:] = a
                    best_cval[si, ri, ti, cvi] = c

                    # finally, predict on the held-out test data here
                    pred = model.predict(tstdat)

                    pred_labs[tstinds] = pred

                if not debug:
                    assert(not np.any(np.isnan(pred_labs)))
                pred_labs_all[ss][ri,tinds] = pred_labs
                
                # compute some performance metrics
                # these are just the "main grid" because easy to compute accuracy
                assert(not np.any(np.isnan(grid_labs_task[is_main_grid_task])))
                acc_bytask[si,ri,ti] = np.mean(pred_labs[is_main_grid_task]==grid_labs_task[is_main_grid_task])
                dprime_bytask[si,ri,ti] = stats_utils.get_dprime(pred_labs[is_main_grid_task], grid_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        save_filename = os.path.join(save_folder, 'decode_multiclass_withintask.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'pred_labs_all': pred_labs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'grid_pts': grid_pts, \
                               'roi_names': roi_names, \
                               })

    
def decode_crosstask(debug=False, n_threads=8, use_bigIPS=False, justIPS=False):
    
    print('debug = %s, n_threads = %d, use_bigIPS = %s, justIPS = %s'%(debug, n_threads, use_bigIPS, justIPS))
    
    if debug:
        subjects=[1]
    else:
        subjects = np.arange(1,8)
    n_subjects = len(subjects)
    n_rois = 11
    make_time_resolved=False

    # first load all data for all subjects, all tasks
    maindat_all = []; repdat_all = []
    mainlabs_all = []; replabs_all = []

    for si, ss in enumerate(subjects):
       
        print('loading S%02d, main task'%ss)
        main_data, _, main_labels, roi_names = data_utils.load_main_task_data(ss, make_time_resolved, use_bigIPS)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            main_data[ri] -= np.tile(np.mean(main_data[ri], axis=1, keepdims=True), [1, main_data[ri].shape[1]])

        maindat_all += [main_data]
        mainlabs_all += [main_labels]
        
        print('loading S%02d, repeat task'%ss)
        rep_data, _, rep_labels, roi_names = data_utils.load_repeat_task_data(ss, make_time_resolved, use_bigIPS)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            rep_data[ri] -= np.tile(np.mean(rep_data[ri], axis=1, keepdims=True), [1, rep_data[ri].shape[1]])

        repdat_all += [rep_data]
        replabs_all += [rep_labels]

    if justIPS:
        n_rois = 1;
        rois_combine = [5,6,7,8]
        print(np.array(roi_names)[rois_combine])
        roi_names = ['IPSall']
        
    # penalties to eval
    c_values = np.logspace(-9, 1, 20)

    n_grid_pts = 16
    
    # store average performance, each roi and subject
    # performance metrics broken down for each task
    n_tasks_trn = 4
    n_tasks_test = 4
    
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks_trn, n_tasks_test))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks_trn, n_tasks_test))
    # performance broken down by points in the grid
    acc_each_point_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_tasks_trn, n_tasks_test)) 
    dprime_each_point_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_tasks_trn, n_tasks_test)) 
    # count individual predictions, to get confusability of categories
    num_preds_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_grid_pts, n_tasks_trn, n_tasks_test)) 
    
    n_cv = 12;
    # store which c value is best
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks_trn, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks_trn, n_cv), np.nan)                        

   
    for si, ss in enumerate(subjects):
        
        main_data = maindat_all[si]
        main_labels = mainlabs_all[si]
        rep_data = repdat_all[si]
        rep_labels = replabs_all[si]

        # gathering labels for main task and for repeat task.
        # all labels will be concatenated [main; repeat]
        inds_use_main = (main_labels['is_main_grid']==True) 
        inds_use_rep = (rep_labels['is_main_grid']==True) 

        xlabs_main = np.array(main_labels['ptx'])[inds_use_main]
        ylabs_main = np.array(main_labels['pty'])[inds_use_main]
        xlabs_rep = np.array(rep_labels['ptx'])[inds_use_rep]
        ylabs_rep = np.array(rep_labels['pty'])[inds_use_rep]
        
        xlabs = np.concatenate([xlabs_main, xlabs_rep], axis=0)
        ylabs = np.concatenate([ylabs_main, ylabs_rep], axis=0)
        
        pt_labs = np.array([xlabs, ylabs]).T
        grid_pts, grid_labs, counts = np.unique(pt_labs, axis=0, return_inverse=True, return_counts=True)
        assert(n_grid_pts==grid_pts.shape[0])

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
        
        is_main_task = task_labs<4
        
        # convert the cv labs, instead of counting over all runs, make them count all runs within 
        # a given task only.
        cv_labs_raw = cv_labs
        cv_labs = np.zeros(np.shape(cv_labs_raw))

        for ti in range(n_tasks_trn):
            inds = task_labs==(ti+1)
            cv_labs[inds] = np.unique(cv_labs_raw[inds], return_inverse=True)[1]

        
        for ri in range(n_rois):
            
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            if justIPS:
                
                # gather data for this ROI
                dat_main = np.concatenate([main_data[ri][inds_use_main,:] \
                           for ri in rois_combine], axis=1)
                dat_rep = np.concatenate([rep_data[ri][inds_use_rep,:] \
                                           for ri in rois_combine], axis=1)
                dat_alltasks = np.concatenate([dat_main, dat_rep], axis=0)
                print(dat_alltasks.shape)
                
            else:

                # gather data for this ROI
                dat_main = main_data[ri][inds_use_main,:]
                dat_rep = rep_data[ri][inds_use_rep,:]
                dat_alltasks = np.concatenate([dat_main, dat_rep], axis=0)

            # loop over tasks to do training on
            for ti, tt in enumerate([1,2,3,4]):
                
                # hold the predicted labels for entire dataset
                pred_labs = np.full(fill_value=np.nan, shape=(np.shape(grid_labs)))

                for cvi, cv in enumerate(np.unique(cv_labs)):

                    if debug & (cvi>0):
                        continue
                        
                    # training dat is the training task of interest, leaving out one run.
                    trninds = (task_labs==tt) & (cv_labs!=cv)

                    trnlabs = grid_labs[trninds]
                    
                    trndat = dat_alltasks[trninds,:]

                    print(' training task %d, cv %d: %d total trials'%(tt, cv, trndat.shape[0]))
                
                    # do regularization parameter (c) selection
                    # this is based on training data only, for the current fold.
                    # cross-validate using leave-one-run-out (for just the training runs here)
                    nest_cv_labs = cv_labs[trninds]
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
                    a = np.mean(model.scores_[0], axis=0)
                    c = model.C_[0]
                    assert(c_values[np.argmax(a)]==c)

                    print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                    sys.stdout.flush()

                    acc_each_cval[si, ri, ti, cvi,:] = a
                    best_cval[si, ri, ti, cvi] = c

                    # loop over tasks to test on
                    for ti2, tt2 in enumerate([1,2,3,4]):
                        
                        # test data is the testing task of interest, left-out run.
                        tstinds = (task_labs==tt2) & (cv_labs==cv)

                        if np.sum(tstinds)==0:
                            print('missing run %d for task %d, skipping'%(cv, tt2))
                            continue
                         
                        tstdat = dat_alltasks[tstinds,:]

                        # finally, predict on the held-out test data here
                        pred = model.predict(tstdat)

                        pred_labs[tstinds] = pred
                        
                if not debug:
                    assert(not np.any(np.isnan(pred_labs)))
                
                # compute some performance metrics
                # these will be over all runs, but within a given task (testing task)
                for ti2, tt2 in enumerate([1,2,3,4]):
                    
                    test_task_inds = (task_labs==tt2)

                    acc_bytask[si,ri,ti,ti2] = np.mean(pred_labs[test_task_inds]==grid_labs[test_task_inds])
                    dprime_bytask[si,ri,ti,ti2] = stats_utils.get_dprime(pred_labs[test_task_inds],
                                                                         grid_labs[test_task_inds])

                    # performance for individual categories (grid pos)
                    for gi, gg in enumerate(np.unique(grid_labs)):

                        # for the trials in categ gg, how often did classifier label correctly (hit rate)?
                        inds = (grid_labs==gg) & test_task_inds 
                        acc_each_point_bytask[si,ri,gi,ti] = np.mean(pred_labs[inds]==grid_labs[inds])
                        
                        # for d-prime, we are going to use all the trials, but re-label them as binary 
                        # for the category of interest. so it captures how well the classifier 
                        # can discriminate this category versus others.
                        pred_labs_binary = (pred_labs==gg).astype(int)
                        grid_labs_binary = (grid_labs==gg).astype(int)
                        dprime_each_point_bytask[si,ri,gi,ti,ti2] = \
                            stats_utils.get_dprime(pred_labs_binary[test_task_inds],
                                                   grid_labs_binary[test_task_inds])

                        # now computing confusability of diff categories 
                        for gi2, gg2 in enumerate(np.unique(grid_labs)):
                            inds = (grid_labs==gg) & (pred_labs==gg2) & test_task_inds
                            num_preds_bytask[si, ri, gi, gi2,ti,ti2] = np.sum(inds)

                        if not debug:
                            assert(np.sum(num_preds_bytask[si,ri,gi,:,ti,ti2])==np.sum((grid_labs==gg) & test_task_inds))

        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        if use_bigIPS:
            save_filename = os.path.join(save_folder, 'decode_multiclass_crosstask_bigIPS.npy')
            if justIPS:
                save_filename = os.path.join(save_folder, 'decode_multiclass_crosstask_IPSall.npy')
        else:
            save_filename = os.path.join(save_folder, 'decode_multiclass_crosstask.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'acc_each_point_bytask': acc_each_point_bytask, \
                               'dprime_each_point_bytask': dprime_each_point_bytask, \
                               'num_preds_bytask': num_preds_bytask, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'grid_pts': grid_pts, \
                               'roi_names': roi_names, \
                               })

        
        

def decode_allmaintask(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    subjects = np.arange(1,8)
    n_subjects = len(subjects)
    n_rois = 11
    make_time_resolved=False

    # first load all data for all subjects
    maindat_all = []; 
    mainlabs_all = []; 

    for si, ss in enumerate(subjects):
       
        print('loading S%02d, main task'%ss)
        main_data, _, main_labels, roi_names = data_utils.load_main_task_data(ss, make_time_resolved)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            main_data[ri] -= np.tile(np.mean(main_data[ri], axis=1, keepdims=True), [1, main_data[ri].shape[1]])

        maindat_all += [main_data]
        mainlabs_all += [main_labels]

    # penalties to eval
    c_values = np.logspace(-9, 1, 20)

    n_grid_pts = 16
    
    # store average performance, each roi and subject
    acc_overall = np.zeros((n_subjects, n_rois))
    dprime_overall = np.zeros((n_subjects, n_rois))
    # performance broken down by points in the grid
    acc_each_point = np.zeros((n_subjects, n_rois, n_grid_pts)) 
    dprime_each_point = np.zeros((n_subjects, n_rois, n_grid_pts)) 

    # count individual predictions, to get confusability of categories
    num_preds = np.zeros((n_subjects, n_rois, n_grid_pts, n_grid_pts)) 
    
    n_cv = 36;
    # store which c value is best
    acc_each_cval = np.full((n_subjects, n_rois, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_cv), np.nan)                        

   
    for si, ss in enumerate(subjects):
        
        main_data = maindat_all[si]
        main_labels = mainlabs_all[si]

        # pull out data from main grid trials only, all tasks here
        inds_use = (main_labels['is_main_grid']==True) 

        # labels are 1-16 for grid positions
        xlabs = np.array(main_labels['ptx'])[inds_use]
        ylabs = np.array(main_labels['pty'])[inds_use]
        pt_labs = np.array([xlabs, ylabs]).T
        grid_pts, grid_labs, counts = np.unique(pt_labs, axis=0, return_inverse=True, return_counts=True)
        assert(n_grid_pts==grid_pts.shape[0])

        # leave-one-run-out
        cv_labs = np.array(main_labels['run_overall'])[inds_use]
        n_cv = len(np.unique(cv_labs))

        for ri in range(n_rois):
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            # data for this ROI
            dat = main_data[ri][inds_use,:]

            # hold the predicted labels for entire dataset
            pred_labs = np.zeros((np.shape(grid_labs)))

            for cvi, cv in enumerate(np.unique(cv_labs)):

                if debug & (cvi>0):
                    continue
                # holding out one run at a time as a test set
                trninds = cv_labs!=cv
                tstinds = cv_labs==cv

                trndat = dat[trninds,:]
                tstdat = dat[tstinds,:]
                trnlabs = grid_labs[trninds]
                tstlabs = grid_labs[tstinds]

                # do regularization parameter (c) selection
                # this is based on training data only, for the current fold.
                # cross-validate using leave-one-run-out (for just the training runs here)
                nest_cv_labs = cv_labs[trninds]
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
                a = np.mean(model.scores_[0], axis=0)
                c = model.C_[0]
                assert(c_values[np.argmax(a)]==c)

                print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                sys.stdout.flush()
                
                acc_each_cval[si, ri, cvi,:] = a
                best_cval[si, ri, cvi] = c

                # finally, predict on the held-out test data here
                pred = model.predict(tstdat)

                pred_labs[tstinds] = pred
                
                
            # compute some performance metrics
            acc_overall[si,ri] = np.mean(pred_labs==grid_labs)
            dprime_overall[si,ri] = stats_utils.get_dprime(pred_labs, grid_labs)

            # performance for individual categories (grid pos)
            for gi, gg in enumerate(np.unique(grid_labs)):

                # for the trials in categ gg, how often did classifier label correctly (hit rate)?
                inds = grid_labs==gg 
                acc_each_point[si,ri,gi] = np.mean(pred_labs[inds]==grid_labs[inds])
                # for d-prime, we are going to use all the trials, but re-label them as binary 
                # for the category of interest. so it captures how well the classifier 
                # can discriminate this category versus others.
                pred_labs_binary = (pred_labs==gg).astype(int)
                grid_labs_binary = (grid_labs==gg).astype(int)
                dprime_each_point[si,ri,gi] = stats_utils.get_dprime(pred_labs_binary, grid_labs_binary)

                # now computing confusability of diff categories 
                for gi2, gg2 in enumerate(np.unique(grid_labs)):
                    inds = (grid_labs==gg) & (pred_labs==gg2)
                    num_preds[si, ri, gi, gi2] = np.sum(inds)

                assert(np.sum(num_preds[si,ri,gi,:])==np.sum(grid_labs==gg))
             
    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    save_filename = os.path.join(save_folder, 'decode_multiclass_allmaintask.npy')
    print('saving to %s'%save_filename)
    np.save(save_filename, {'acc_overall': acc_overall, \
                           'dprime_overall': dprime_overall, \
                           'acc_each_point': acc_each_point, \
                           'dprime_each_point': dprime_each_point, \
                           'num_preds': num_preds, \
                           'acc_each_cval': acc_each_cval, \
                           'best_cval': best_cval, \
                           })


def decode_alltasks(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    subjects = np.arange(1,8)
    n_subjects = len(subjects)
    n_rois = 11
    make_time_resolved=False

    # first load all data for all subjects, all tasks
    maindat_all = []; repdat_all = []
    mainlabs_all = []; replabs_all = []

    for si, ss in enumerate(subjects):
       
        print('loading S%02d, main task'%ss)
        main_data, _, main_labels, roi_names = data_utils.load_main_task_data(ss, make_time_resolved)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            main_data[ri] -= np.tile(np.mean(main_data[ri], axis=1, keepdims=True), [1, main_data[ri].shape[1]])

        maindat_all += [main_data]
        mainlabs_all += [main_labels]
        
        print('loading S%02d, repeat task'%ss)
        rep_data, _, rep_labels, roi_names = data_utils.load_repeat_task_data(ss, make_time_resolved)

        for ri in range(n_rois):
            # subtract mean across voxels each trial
            rep_data[ri] -= np.tile(np.mean(rep_data[ri], axis=1, keepdims=True), [1, rep_data[ri].shape[1]])

        repdat_all += [rep_data]
        replabs_all += [rep_labels]

    # penalties to eval
    c_values = np.logspace(-9, 1, 20)

    n_grid_pts = 16
    
    # store average performance, each roi and subject
    acc_overall = np.zeros((n_subjects, n_rois))
    dprime_overall = np.zeros((n_subjects, n_rois))
    # performance broken down by points in the grid
    acc_each_point = np.zeros((n_subjects, n_rois, n_grid_pts)) 
    dprime_each_point = np.zeros((n_subjects, n_rois, n_grid_pts)) 
    # count individual predictions, to get confusability of categories
    num_preds = np.zeros((n_subjects, n_rois, n_grid_pts, n_grid_pts)) 
    
    # performance metrics broken down for each task
    n_tasks = 4
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks))
    # performance broken down by points in the grid
    acc_each_point_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_tasks)) 
    dprime_each_point_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_tasks)) 
    # count individual predictions, to get confusability of categories
    num_preds_bytask = np.zeros((n_subjects, n_rois, n_grid_pts, n_grid_pts, n_tasks)) 
    
    n_cv = 48;
    # store which c value is best
    acc_each_cval = np.full((n_subjects, n_rois, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_cv), np.nan)                        

   
    for si, ss in enumerate(subjects):
        
        main_data = maindat_all[si]
        main_labels = mainlabs_all[si]
        rep_data = repdat_all[si]
        rep_labels = replabs_all[si]

        # gathering labels for main task and for repeat task.
        # all labels will be concatenated [main; repeat]
        inds_use_main = (main_labels['is_main_grid']==True) 
        inds_use_rep = (rep_labels['is_main_grid']==True) 

        xlabs_main = np.array(main_labels['ptx'])[inds_use_main]
        ylabs_main = np.array(main_labels['pty'])[inds_use_main]
        xlabs_rep = np.array(rep_labels['ptx'])[inds_use_rep]
        ylabs_rep = np.array(rep_labels['pty'])[inds_use_rep]
        
        xlabs = np.concatenate([xlabs_main, xlabs_rep], axis=0)
        ylabs = np.concatenate([ylabs_main, ylabs_rep], axis=0)
        
        pt_labs = np.array([xlabs, ylabs]).T
        grid_pts, grid_labs, counts = np.unique(pt_labs, axis=0, return_inverse=True, return_counts=True)
        assert(n_grid_pts==grid_pts.shape[0])

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
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            # data for this ROI
            dat_main = main_data[ri][inds_use_main,:]
            dat_rep = rep_data[ri][inds_use_rep,:]
            dat = np.concatenate([dat_main, dat_rep], axis=0)
            
            # hold the predicted labels for entire dataset
            pred_labs = np.zeros((np.shape(grid_labs)))

            for cvi, cv in enumerate(np.unique(cv_labs)):

                if debug & (cvi>0):
                    continue
                # holding out one run at a time as a test set
                trninds = cv_labs!=cv
                tstinds = cv_labs==cv

                trndat = dat[trninds,:]
                tstdat = dat[tstinds,:]
                trnlabs = grid_labs[trninds]
                tstlabs = grid_labs[tstinds]

                # do regularization parameter (c) selection
                # this is based on training data only, for the current fold.
                # cross-validate using leave-one-run-out (for just the training runs here)
                nest_cv_labs = cv_labs[trninds]
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
                a = np.mean(model.scores_[0], axis=0)
                c = model.C_[0]
                assert(c_values[np.argmax(a)]==c)

                print('    cv fold %d (elapsed = %.6f s): best c = %.5f, max acc = %.2f'%(cvi, elapsed, c, np.max(a)))
                sys.stdout.flush()
                
                acc_each_cval[si, ri, cvi,:] = a
                best_cval[si, ri, cvi] = c

                # finally, predict on the held-out test data here
                pred = model.predict(tstdat)

                pred_labs[tstinds] = pred
                
                
            # compute some performance metrics
            acc_overall[si,ri] = np.mean(pred_labs==grid_labs)
            dprime_overall[si,ri] = stats_utils.get_dprime(pred_labs, grid_labs)

            # performance for individual categories (grid pos)
            for gi, gg in enumerate(np.unique(grid_labs)):

                # for the trials in categ gg, how often did classifier label correctly (hit rate)?
                inds = grid_labs==gg 
                acc_each_point[si,ri,gi] = np.mean(pred_labs[inds]==grid_labs[inds])
                # for d-prime, we are going to use all the trials, but re-label them as binary 
                # for the category of interest. so it captures how well the classifier 
                # can discriminate this category versus others.
                pred_labs_binary = (pred_labs==gg).astype(int)
                grid_labs_binary = (grid_labs==gg).astype(int)
                dprime_each_point[si,ri,gi] = stats_utils.get_dprime(pred_labs_binary, grid_labs_binary)

                # now computing confusability of diff categories 
                for gi2, gg2 in enumerate(np.unique(grid_labs)):
                    inds = (grid_labs==gg) & (pred_labs==gg2)
                    num_preds[si, ri, gi, gi2] = np.sum(inds)

                assert(np.sum(num_preds[si,ri,gi,:])==np.sum(grid_labs==gg))
             
            # now looking at performance for each task individually
            for ti, tt in enumerate([1,2,3,4]):

                tinds = task_labs==tt
                print(np.sum(tinds))
                
                acc_bytask[si,ri,ti] = np.mean(pred_labs[tinds]==grid_labs[tinds])
                dprime_bytask[si,ri,ti] = stats_utils.get_dprime(pred_labs[tinds], grid_labs[tinds])

                # performance for individual categories (grid pos)
                for gi, gg in enumerate(np.unique(grid_labs)):

                    # for the trials in categ gg, how often did classifier label correctly (hit rate)?
                    inds = tinds & (grid_labs==gg)
                    acc_each_point_bytask[si,ri,gi,ti] = np.mean(pred_labs[inds]==grid_labs[inds])
                    # for d-prime, we are going to use all the trials, but re-label them as binary 
                    # for the category of interest. so it captures how well the classifier 
                    # can discriminate this category versus others.
                    pred_labs_binary = (pred_labs[tinds]==gg).astype(int)
                    grid_labs_binary = (grid_labs[tinds]==gg).astype(int)
                    dprime_each_point_bytask[si,ri,gi,ti] = stats_utils.get_dprime(pred_labs_binary, grid_labs_binary)

                    # now computing confusability of diff categories 
                    for gi2, gg2 in enumerate(np.unique(grid_labs)):
                        inds = (grid_labs[tinds]==gg) & (pred_labs[tinds]==gg2)
                        num_preds_bytask[si, ri, gi, gi2, ti] = np.sum(inds)

                    assert(np.sum(num_preds_bytask[si,ri,gi,:,ti])==np.sum(grid_labs[tinds]==gg))
             
    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    save_filename = os.path.join(save_folder, 'decode_multiclass_alltasks.npy')
    print('saving to %s'%save_filename)
    np.save(save_filename, {'acc_overall': acc_overall, \
                           'dprime_overall': dprime_overall, \
                           'acc_each_point': acc_each_point, \
                           'dprime_each_point': dprime_each_point, \
                           'num_preds': num_preds, \
                           'acc_bytask': acc_bytask, \
                           'dprime_bytask': dprime_bytask, \
                           'acc_each_point_bytask': acc_each_point_bytask, \
                           'dprime_each_point_bytask': dprime_each_point_bytask, \
                           'num_preds_bytask': num_preds_bytask, \
                           'acc_each_cval': acc_each_cval, \
                           'best_cval': best_cval, \
                           'grid_pts': grid_pts, \
                           'roi_names': roi_names, \
                           })


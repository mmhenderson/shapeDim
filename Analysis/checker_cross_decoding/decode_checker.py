import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis, sklearn.linear_model
import time


# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import decoding_utils
from code_utils import stats_utils, grid_utils


def decode_withintask(debug=False, n_threads=8):
    
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
    # penalties to eval
    c_values = np.logspace(-9, 1, 20)
    
    n_tasks = 4
    
    n_trnsides = 4 # which half of the grid to train on?
    trnsides = ['left','right','bottom','top']
    
    # store average performance, each roi and subject
    # these will all be for just "main grid" trials
    acc_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_trnsides))
    dprime_bytask = np.zeros((n_subjects, n_rois, n_tasks, n_trnsides))

    n_cv = 12;
    # store which c value is best
    # these accs are on the nested training data, not the testing data
    acc_each_cval = np.full((n_subjects, n_rois, n_tasks, n_trnsides, n_cv, len(c_values)), np.nan)
    best_cval = np.full((n_subjects, n_rois, n_tasks, n_trnsides, n_cv), np.nan)                        
   
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
        n_cv = len(np.unique(cv_labs))
        
        pt_labs = np.array([concat_labels['ptx'], concat_labels['pty']]).T
        
        # category is always checkerboard task here
        categ_labs = grid_utils.get_categ(pt_labs, task_num=3)
        
        is_main_grid = np.array(concat_labels['is_main_grid']==True) 

    
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
                
                categ_labs_task = categ_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                is_main_grid_task = is_main_grid[tinds]
                pt_labs_task = pt_labs[tinds,:]
                
                # data for this ROI
                data_task = data[tinds,:]
                dat = data_task
                
                print(' processing task %d, %d total trials'%(tt, dat.shape[0]))

                for tri in range(n_trnsides):
                    
                    # train on one half of the shape space
                    # so we can test generalization to other side
                    if trnsides[tri]=='left':
                        pt_inds_trn = pt_labs_task[:,0]<2.5
                    elif trnsides[tri]=='right':
                        pt_inds_trn = pt_labs_task[:,0]>2.5
                    elif trnsides[tri]=='bottom':
                        pt_inds_trn = pt_labs_task[:,1]<2.5
                    elif trnsides[tri]=='top':
                        pt_inds_trn = pt_labs_task[:,1]>2.5
                        
                    print('   %s, %d trials'%(trnsides[tri], np.sum(pt_inds_trn)))
                    
                    # hold the predicted labels for this task
                    nt = len(categ_labs_task)
                    pred_labs = np.full(fill_value=np.nan, shape=[nt,])
                    prob_each = np.full(fill_value=np.nan, shape=[nt,2])

                    for cvi, cv in enumerate(np.unique(cv_labs_task)):

                        if debug & (cvi>0):
                            continue

                        # holding out one run at a time as a test set
                        # training set is all the other runs, only main grid trials.
                        trninds = (cv_labs_task!=cv) & (is_main_grid_task) & (pt_inds_trn)
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

                        acc_each_cval[si, ri, ti, tri, cvi,:] = a
                        best_cval[si, ri, ti, tri, cvi] = c

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
                    preds_all[si][ri][ti][tri] = pred_labs
                    probs_all[si][ri][ti][tri] = prob_each

                    # compute some performance metrics
                    # these are just the "main grid" because easy to compute accuracy
                    assert(not np.any(np.isnan(categ_labs_task[is_main_grid_task])))
                    acc_bytask[si,ri,ti,tri] = \
                        np.mean(pred_labs[is_main_grid_task]==categ_labs_task[is_main_grid_task])
                    dprime_bytask[si,ri,ti,tri] = \
                        stats_utils.get_dprime(pred_labs[is_main_grid_task], categ_labs_task[is_main_grid_task])

        # save after each subject, in case of a crash
        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        save_filename = os.path.join(save_folder, 'decode_checker_across_sides.npy')
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_bytask': acc_bytask, \
                               'dprime_bytask': dprime_bytask, \
                               'preds_all': preds_all, \
                               'probs_all': probs_all, \
                               'acc_each_cval': acc_each_cval, \
                               'best_cval': best_cval, \
                               'roi_names': roi_names, \
                               })

        
import numpy as np
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
from code_utils import stats_utils

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

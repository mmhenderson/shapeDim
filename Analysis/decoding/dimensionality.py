import numpy as np
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis, sklearn.linear_model
import time
import itertools

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])
# root = /usr/local/serenceslab/maggie/shapeDim/

from code_utils import file_utils, data_utils
from code_utils import decoding_utils
from code_utils import stats_utils

def compute_shattering_dim(debug=False, n_threads=8):
    
    print('debug = %s, n_threads = %d'%(debug, n_threads))
    
    # subjects=[1,2,3,4]
    subjects=[5,6,7]
    # subjects = np.arange(1,8)
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
    c = 0.005
    
    n_grid_pts = 16
    
    # define here all the binary dichotomies that we will be testing
    dich = get_all_dichotomies(n_grid_pts)
    n_dich = dich.shape[0] # for 16, should be 6435
    
    # store the decoding performance, for each dichotomy
    n_tasks = 4
    
    for si, ss in enumerate(subjects):
        
        acc_each_dichotomy = np.zeros((n_rois, n_tasks, n_dich))
    
        
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
        
        for ri in range(n_rois):
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            for ti, tt in enumerate([1,2,3,4]):
                
                tinds = task_labs==tt
                
                grid_labs_task = grid_labs[tinds]
                cv_labs_task = cv_labs[tinds]
                
                # data for this ROI
                dat_main = main_data[ri][inds_use_main,:][tinds[is_main_task],:]
                dat_rep = rep_data[ri][inds_use_rep,:][tinds[~is_main_task],:]
                dat = np.concatenate([dat_main, dat_rep], axis=0)

                print('processing task %d: %d total trials'%(tt, dat.shape[0]))
                
                for di in range(n_dich):
    
                    if debug and (di>1):
                        continue
                
                    # create binarized labels for this dichotomy
                    d = dich[di,:]
                    pts1 = np.arange(n_grid_pts)[d==1]
                    print('grouping:')
                    print(pts1)
                    labels_dich = np.isin(grid_labs_task, pts1).astype(int)
                    assert(np.mean(labels_dich)==1/2)

                    st = time.time()

                    # do the decoding here
                    # set up cross-validation, leave-one-run-out
                    # using this generator is faster than manually looping
                    cv_obj = sklearn.model_selection.LeaveOneGroupOut()
                    cv_generator = cv_obj.split(dat, labels_dich, cv_labs_task)
                    # define model
                    model = sklearn.linear_model.LogisticRegressionCV(Cs = [c], \
                                                                      cv = cv_generator, \
                                                                    solver='lbfgs', \
                                                                    penalty='l2', \
                                                                    n_jobs = n_threads , \
                                                                    max_iter = 1000)
                    model.fit(dat, labels_dich)
                    acc = np.mean(model.scores_[1])

                    elapsed = time.time() - st
                    print('acc = %.2f, elapsed = %.5f s'%(acc, elapsed))

                    acc_each_dichotomy[ri, ti, di] = acc
    

        save_folder = os.path.join(root, 'Analysis', 'decoding_results')
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
        # save_filename = os.path.join(save_folder, 'shattering_dim.npy')
        save_filename = os.path.join(save_folder, 'S%02d_shattering_dim.npy'%ss)
        print('saving to %s'%save_filename)
        np.save(save_filename, {'acc_each_dichotomy': acc_each_dichotomy, \
                                'c_value': c, \
                               'grid_pts': grid_pts, \
                               'roi_names': roi_names, \
                               })    
    
def get_all_dichotomies(n_pts):
    
    assert(np.mod(n_pts,2)==0) # has to be even
    
    # list all the possible ways of taking half the points
    n_half = int(n_pts/2)
    groups = list(itertools.combinations(np.arange(n_pts),n_half))
    dich = []
    for gi, gg in enumerate(groups):
        d = np.isin(np.arange(n_pts), gg).astype(int)
        d = d[None,:]
        if gi==0:
            dich = d
        else:
            # check if this split or its inverse is already in my list
            already_in = np.any(np.all(d==dich, axis=1)) or (np.any(np.all((1-d)==dich, axis=1)))
            if not already_in:
                dich = np.concatenate([dich, d], axis=0)
    
    return dich


def compute_ccgp(debug=False, n_threads=8):
    
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
    c = 0.005
    

    # store the decoding performance, for each dichotomy
    n_tasks = 4
    acc_ccgp = np.zeros((n_subjects, n_rois, n_tasks, 4))
    
    for si, ss in enumerate(subjects):
        
        main_data = maindat_all[si]
        main_labels = mainlabs_all[si]
        rep_data = repdat_all[si]
        rep_labels = replabs_all[si]

        # gathering labels for main task and for repeat task.
        # all labels will be concatenated [main; repeat]
        inds_use_main = (main_labels['is_main_grid']==True) 
        inds_use_rep = (rep_labels['is_main_grid']==True) 

        center = 2.5
        # binarize these, two categories on each axis
        xlabs_main = (np.array(main_labels['ptx'])[inds_use_main]>center).astype(int)
        ylabs_main = (np.array(main_labels['pty'])[inds_use_main]>center).astype(int)
        xlabs_rep = (np.array(rep_labels['ptx'])[inds_use_rep]>center).astype(int)
        ylabs_rep = (np.array(rep_labels['pty'])[inds_use_rep]>center).astype(int)
        
        xlabs = np.concatenate([xlabs_main, xlabs_rep], axis=0)
        ylabs = np.concatenate([ylabs_main, ylabs_rep], axis=0)
        
        # repeat task is task "4" out of 4 here
        task_labs_main = np.array(main_labels['task'])[inds_use_main]
        task_labs_rep = 4 * np.ones((np.sum(inds_use_rep), ), dtype=int)
        task_labs = np.concatenate([task_labs_main, task_labs_rep], axis=0)

        is_main_task = task_labs<4
        
        for ri in range(n_rois):
            
            if debug & (ri>0):
                continue
            print('proc S%02d, %s'%(ss, roi_names[ri]))
            
            for ti, tt in enumerate([1,2,3,4]):
                
                tinds = task_labs==tt
                
                xlabs_task = xlabs[tinds]
                ylabs_task = ylabs[tinds]
               
                # data for this ROI
                dat_main = main_data[ri][inds_use_main,:][tinds[is_main_task],:]
                dat_rep = rep_data[ri][inds_use_rep,:][tinds[~is_main_task],:]
                dat = np.concatenate([dat_main, dat_rep], axis=0)

                print('processing task %d: %d total trials'%(tt, dat.shape[0]))
                
                xx = -1
                # splitting by x or y, decoding other dimension
                for split_dim, dec_dim in zip([xlabs_task, ylabs_task], [ylabs_task, xlabs_task]):
                    
                    # loop over training/testing splits
                    for trni, testi in zip([0,1], [1,0]):
                        
                        # split the data according to this dimension
                        trninds = split_dim==trni
                        tstinds = split_dim==testi
                        
                        # define model
                        st = time.time()
                        model = sklearn.linear_model.LogisticRegression(C = c, \
                                                                        solver='lbfgs', \
                                                                        penalty='l2', \
                                                                        n_jobs = n_threads , \
                                                                        max_iter = 1000)
                        model.fit(dat[trninds,:], dec_dim[trninds])
                        acc = model.score(dat[tstinds,:], dec_dim[tstinds])
                        
                        elapsed = time.time() - st
                        print('acc = %.2f, elapsed = %.5f s'%(acc, elapsed))

                        xx+=1
                        acc_ccgp[si,ri,ti,xx] = acc

    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    save_filename = os.path.join(save_folder, 'ccgp.npy')
    print('saving to %s'%save_filename)
    np.save(save_filename, {'acc_ccgp': acc_ccgp, \
                            'c_value': c, \
                           'roi_names': roi_names, \
                           })

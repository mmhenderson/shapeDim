import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os, sys
import pandas as pd
import sklearn
import sklearn.svm, sklearn.discriminant_analysis
import time

root = '/usr/local/serenceslab/maggie/shapeDim/'

sys.path.append(os.path.join(root, 'Analysis'))
from code_utils import file_utils, data_utils
from code_utils import decoding_utils

def decode(subjects = np.arange(1,8)):

    n_subj = len(subjects)

    task_names = ['Linear (1)','Linear (2)','Checker'];
    n_tasks = len(task_names)

    # three different ways to do binary decoding
    n_bounds = 3;
    bound_names = ['Decode: Linear (1)','Decode: Linear (2)','Decode: Checker'];
    quad_groups = [[[1, 4], [2, 3]],
                    [[1, 2], [3, 4]],
                    [[1, 3], [2, 4]]];

    make_time_resolved = False

    n_rois = 11
    
    # first load all data for all subjects, both tasks
    maindat_all = []; repdat_all = []
    mainlabs_all = []; replabs_all = []

    for si, ss in enumerate(subjects):
        # si = 0; ss = 1;

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

    # now do the decoding, one subject/ROI at a time
    
    # decode_func = decode_func_svc
    decode_func = decoding_utils.decode_func_lda
    # decode_func = decode_func_normeucdist

    acc = np.zeros((n_subj, n_rois, n_tasks+1, n_bounds))

    for si, ss in enumerate(subjects):

        main_data = maindat_all[si]
        main_labels = mainlabs_all[si]

        rep_data = repdat_all[si]
        rep_labels = replabs_all[si]

        for ri in range(n_rois):

            print('proc S%02d, %s'%(ss, roi_names[ri]))
            st = time.time()

            for ti in range(n_tasks):

                # pull out data from main grid trials only, this task
                inds_use = (main_labels['is_main_grid']==True) & (main_labels['task']==ti+1)

                dat = main_data[ri][inds_use,:]
                quadlabs = np.array(main_labels['quadrant'])[inds_use]

                cv_labs = np.array(main_labels['run_overall'])[inds_use]
                n_cv = len(np.unique(cv_labs))

                for bi in range(n_bounds):

                    # convert quadrants to binary labels, for the current boundary
                    categ_labs = np.zeros(np.shape(quadlabs))
                    categ_labs[np.isin(quadlabs, quad_groups[bi][0])] = 1
                    categ_labs[np.isin(quadlabs, quad_groups[bi][1])] = 2
                    assert(not np.any(categ_labs==0))

                    pred_categ_labs = np.zeros(np.shape(quadlabs))

                    for cvi, cv in enumerate(np.unique(cv_labs)):

                        trninds = cv_labs!=cv
                        tstinds = cv_labs==cv

                        # check for balanced labels
                        un, counts = np.unique(categ_labs[trninds], return_counts=True)
                        assert(counts[0]==counts[1])

                        # st = time.time()
                        pred = decode_func(dat[trninds,:], categ_labs[trninds], dat[tstinds,:])
                        # elapsed = time.time() - st
                        # print('elapsed time: %.5f s'%elapsed)

                        pred_categ_labs[tstinds] = pred

                    acc[si, ri, ti, bi] = np.mean(categ_labs==pred_categ_labs)


            # repeat task is the "4th" task here
            ti = 3;

            # pull out data from main grid trials only, this task
            inds_use = (rep_labels['is_main_grid']==True)
            dat = rep_data[ri][inds_use,:]
            quadlabs = np.array(rep_labels['quadrant'])[inds_use]

            cv_labs = np.array(rep_labels['run_overall'])[inds_use]
            n_cv = len(np.unique(cv_labs))

            for bi in range(n_bounds):

                # convert quadrants to binary labels, for the current boundary
                categ_labs = np.zeros(np.shape(quadlabs))
                categ_labs[np.isin(quadlabs, quad_groups[bi][0])] = 1
                categ_labs[np.isin(quadlabs, quad_groups[bi][1])] = 2
                assert(not np.any(categ_labs==0))

                pred_categ_labs = np.zeros(np.shape(quadlabs))

                for cvi, cv in enumerate(np.unique(cv_labs)):

                    trninds = cv_labs!=cv
                    tstinds = cv_labs==cv

                    # check for balanced labels
                    un, counts = np.unique(categ_labs[trninds], return_counts=True)
                    assert(counts[0]==counts[1])

                    # st = time.time()
                    pred = decode_func(dat[trninds,:], categ_labs[trninds], dat[tstinds,:])
                    # elapsed = time.time() - st
                    # print('elapsed time: %.5f s'%elapsed)

                    pred_categ_labs[tstinds] = pred

                acc[si, ri, ti, bi] = np.mean(categ_labs==pred_categ_labs)

            elapsed = time.time() - st
            print('elapsed time: %.5f s'%elapsed)
            
            
    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    save_filename = os.path.join(save_folder, 'decode_binary_within_task.npy')
    np.save(save_filename, acc)

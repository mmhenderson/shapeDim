import os, sys
import numpy as np
import pandas as pd

path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])

sys.path.append(os.path.join(root, 'Analysis'))
from code_utils import stats_utils, grid_utils, data_utils


def bootstrap_binary_hardtrials(n_boot_iter = 1000, n_boot_samp = 50, \
                                rndseed = 243555, correct_only = True):

    # computing the classifier confidence for "hard" trials in each binary task
    # based on the predictions of multinomial classifier (from decode_multiclass.py)
    # using bootstrap resampling to equate trial numbers/coordinates across tasks
    
    # load results of multinomial classifier
    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    save_filename = os.path.join(save_folder, 'decode_multiclass_withintask.npy')
    dec_withintask = np.load(save_filename, allow_pickle=True).item()

    # specify some parameters
    roi_names = dec_withintask['roi_names']
    n_rois = len(roi_names)
    n_tasks = 4;
    n_subjects = 7
    subjects = np.arange(1,8)
    
    # load labels for each trial
    lab = dict()
    for ss in subjects:

        # get labels for all the trials, this subject
        main_labels = data_utils.load_main_task_labels(ss)
        rep_labels = data_utils.load_repeat_task_labels(ss)
        lab[ss] = pd.concat([main_labels, rep_labels], axis=0)

        
    grid_pts = grid_utils.get_main_grid()
    # NOTE i am swapping the columns here
    # because this is the order you get from doing np.unique(pts)
    # this is the actual order that the predictions 1-16 of this classifier
    # correspond to. 
    grid_pts = grid_pts[:,[1,0]] 
    

    # create a set of bins for the shape-space coordinates that are sampled
    # these only span the "hard" part of coordinate space
    n_coord_bins = 6;
    # n_coord_bins = 4;
    coord_bin_edges = np.linspace(2.5-0.701, 2.5+0.701, n_coord_bins+1)
    center=2.5
    bin_centers = coord_bin_edges[0:-1]+(coord_bin_edges[1]-coord_bin_edges[0])/2
    bin_dist = (center-bin_centers).round(2)

    
    # [subjects, rois, tasks, axes, bootstrap iterations]
    signedconf_hardtrials_boot = np.zeros((n_subjects, n_rois, 2, 2, n_boot_iter))
    dprime_hardtrials_boot = np.zeros((n_subjects, n_rois, 2, 2, n_boot_iter))

    np.random.seed(rndseed)

    for si, ss in enumerate(subjects):

        print(si)

        # loop over "axes" - which boundary to compute accuracy for
        for ii in [0,1]:

            l = lab[ss]

            task_labs = l['task']

            pt_labs = np.array([l['ptx'], l['pty']]).T

            is_main_grid = l['is_main_grid']==1

            correct = np.array(l['subject_correct'])

            # using only the trials in center of grid
            # (these are hard for all the tasks)
            dist_from_center1 = l['dist_from_bound1']
            dist_from_center2 = l['dist_from_bound2']
            is_hard = (dist_from_center1<0.79) & (dist_from_center2<0.79) & (~is_main_grid)

            # actual label along the axis of interest
            categ_actual = np.array(l['categ_task%d'%(ii+1)]).astype(int)

            # actual coordinate along the axis of interest
            coord_actual = pt_labs[:,ii].round(2)

            # break the coordinate axis into bins, use these to balance the stimulus
            # distribution when bootstrapping
            coord_binned = np.zeros(coord_actual.shape, dtype=int)-1
            for cb in np.arange(n_coord_bins):
                binds = (coord_actual>=coord_bin_edges[cb]) & (coord_actual<=coord_bin_edges[cb+1])
                coord_binned[binds] = cb
            assert(np.all(coord_binned[is_hard]>-1))

            if correct_only:
                inds1 = np.where(is_hard & correct & (task_labs==1))[0]
                inds2 = np.where(is_hard & correct & (task_labs==2))[0]
            else:
                inds1 = np.where(is_hard & (task_labs==1))[0]
                inds2 = np.where(is_hard & (task_labs==2))[0]

            
            # figure out which bins we can use and still have everything balanced in both tasks
            un1, counts1 = np.unique(coord_binned[inds1], return_counts=True)
            un2, counts2 = np.unique(coord_binned[inds2], return_counts=True)

            # print(un1, counts1)
            # print(un2, counts2)
            # print(bin_dist[un1], bin_dist[un2])

            bins_balance = []
            for uu in np.union1d(un1, un2):
                d = bin_dist[uu]
                # to use this this bin, both this bin and the one opposite it (across the 
                # boundary) have to be sampled in both tasks
                in1 = (d in bin_dist[un1]) and (-d in bin_dist[un1])
                in2 = (d in bin_dist[un2]) and (-d in bin_dist[un2])
                if in1 and in2:
                    bins_balance += [uu]

            # print(bin_dist[bins_balance])

            # checking that the bins we are using represent each category equally
            assert(np.mean(bin_dist[bins_balance]<0)==0.5)

            # decide how many samples per bin, based on how many bins we have
            n_samp_eachbin = int(np.ceil(n_boot_samp/len(bins_balance)))

            is_hard_all = is_hard
            coord_binned_all = coord_binned
            categ_actual_all = categ_actual
            correct_all = correct

            # loop over the tasks
            for ti,tt in enumerate([1,2]):

                is_hard = is_hard_all[task_labs==tt]
                coord_binned = coord_binned_all[task_labs==tt]
                categ_actual = categ_actual_all[task_labs==tt]
                correct = correct_all[task_labs==tt]

                if correct_only:
                    inds = np.where(is_hard & correct)[0]
                else:
                    inds = np.where(is_hard)[0]
                    
                for bi in range(n_boot_iter):

                    # make a resampling order that represents each bin equally
                    inds_resamp = []
                    for bn in bins_balance:
                        inds_bin = inds[coord_binned[inds]==bn]
                        assert(len(inds_bin)>0)
                        # if bi==0:
                        #     print(len(inds_bin), n_samp_eachbin)
                        inds_resamp.append(np.random.choice(inds_bin, n_samp_eachbin, replace=True))    
                    inds_resamp = np.concatenate(inds_resamp, axis=0)

                    # double check resample order
                    assert(np.mean(categ_actual[inds_resamp]==1)==0.5)
                    assert(np.all(np.isin(coord_binned[inds_resamp], bins_balance)))
                    counts = np.array([np.sum(coord_binned[inds_resamp]==bn) for bn in bins_balance])
                    assert(np.all(counts==n_samp_eachbin))

                    # loop over ROIs
                    for ri in range(n_rois):

                        pred = dec_withintask['preds_all'][si][ri][ti].astype(int)

                        # which binary category did the classifier predict?
                        categ_pred = grid_utils.get_categ(grid_pts[pred,:], ii+1)

                        prob = dec_withintask['probs_all'][si][ri][ti]

                        # "confidence" in assignment to category 2 vs 1
                        # group the 16 points into categories w/r/t relevant axis
                        g1 = grid_utils.get_categ(grid_pts, ii+1)==1
                        p_categ1 = np.sum(prob[:,g1], axis=1)
                        g2 = grid_utils.get_categ(grid_pts, ii+1)==2
                        p_categ2 = np.sum(prob[:,g2], axis=1)

                        # signed confidence will be: p(correct) - p(incorrect)
                        signedconf = np.zeros_like(p_categ1)
                        signedconf[categ_actual==1] = p_categ1[categ_actual==1] - p_categ2[categ_actual==1]
                        signedconf[categ_actual==2] = p_categ2[categ_actual==2] - p_categ1[categ_actual==2]                    


                        d = stats_utils.get_dprime(categ_pred[inds_resamp], categ_actual[inds_resamp])
                        dprime_hardtrials_boot[si,ri,ti,ii,bi] = d;

                        signedconf_hardtrials_boot[si,ri,ti,ii,bi] = np.mean(signedconf[inds_resamp])
        

    # save this, just because it takes a long time to run
    if correct_only:
        fn2save = os.path.join(save_folder, 'decode_multiclass_binary_hardtrials_bootstrap_correctonly.npy')
    else:
        fn2save = os.path.join(save_folder, 'decode_multiclass_binary_hardtrials_bootstrap.npy')
        
    np.save(fn2save, {'signedconf_hardtrials_boot': signedconf_hardtrials_boot, \
                      'dprime_hardtrials_boot': dprime_hardtrials_boot, \
                     })
    
    
    
    
def bootstrap_correct_incorrect(n_boot_iter = 1000, n_boot_samp = 50, rndseed = 546466):
    
    
    # computing the classifier confidence for correct vs. incorrect trials each task
    # based on the predictions of multinomial classifier (from decode_multiclass.py)
    # using bootstrap resampling to equate trial numbers/coordinates across tasks
    
    # load results of multinomial classifier
    save_folder = os.path.join(root, 'Analysis', 'decoding_results')
    save_filename = os.path.join(save_folder, 'decode_multiclass_withintask.npy')
    dec_withintask = np.load(save_filename, allow_pickle=True).item()

    # specify some parameters
    roi_names = dec_withintask['roi_names']
    n_rois = len(roi_names)
    n_tasks = 4;
    n_subjects = 7
    subjects = np.arange(1,8)
    
    # load labels for each trial
    lab = dict()
    for ss in subjects:

        # get labels for all the trials, this subject
        main_labels = data_utils.load_main_task_labels(ss)
        rep_labels = data_utils.load_repeat_task_labels(ss)
        lab[ss] = pd.concat([main_labels, rep_labels], axis=0)

        
    grid_pts = grid_utils.get_main_grid()
    # NOTE i am swapping the columns here
    # because this is the order you get from doing np.unique(pts)
    # this is the actual order that the predictions 1-16 of this classifier
    # correspond to. 
    grid_pts = grid_pts[:,[1,0]] 
    
    
    
    n_coord_bins = 12;
    coord_bin_edges = np.linspace(-0.701, 0.701, n_coord_bins+1)
    center=2.5
    bin_centers = coord_bin_edges[0:-1]+(coord_bin_edges[1]-coord_bin_edges[0])/2
    bin_dist = bin_centers.round(2)

    # [subjects, rois, tasks, correct/incorrect, bootstrap iterations]
    signedconf_hardtrials_sepcorrect_boot = np.zeros((n_subjects, n_rois, 3, 2, n_boot_iter))
    dprime_hardtrials_sepcorrect_boot = np.zeros((n_subjects, n_rois, 3, 2, n_boot_iter))

    np.random.seed(rndseed)

    for si, ss in enumerate(subjects):

        print(si)

        for ti, tt in enumerate([1,2,3]):

            l = lab[ss][lab[ss]['task']==tt]

            pt_labs = np.array([l['ptx'], l['pty']]).T

            is_main_grid = l['is_main_grid']==1

            ii = ti; # focusing on the task-relevant axis here

            # is it a hard trial?
            is_hard = ~is_main_grid


            categ_actual = np.array(l['categ_task%d'%(ii+1)])

            coord_actual = np.array(l['dist_from_bound%d'%(ii+1)])
            coord_actual[categ_actual==1] = (-1)*coord_actual[categ_actual==1]

            # coord_actual = pt_labs[:,ii].round(2)

            # make rough bins that help equate stimulus distributions
            # across tasks when bootstrapping
            coord_binned = np.zeros(coord_actual.shape, dtype=int)-1
            for cb in np.arange(n_coord_bins):
                binds = (coord_actual>=coord_bin_edges[cb]) & (coord_actual<=coord_bin_edges[cb+1])
                coord_binned[binds] = cb
            assert(np.all(coord_binned[is_hard]>-1))


            # was the subject correct or incorrect?
            correct = np.array(l['subject_correct'])

            inds1 = np.where(is_hard & correct)[0]
            inds2 = np.where(is_hard & ~correct)[0]

            # print(len(inds1), len(inds2))

            # now figure out which bins we can use and still have everything balanced in both correct/incorrect
            un1, counts1 = np.unique(coord_binned[inds1], return_counts=True)
            un2, counts2 = np.unique(coord_binned[inds2], return_counts=True)

            # print(un1, counts1)
            # print(un2, counts2)

            # print(bin_dist[un1], bin_dist[un2])

            bins_balance = []
            for uu in np.union1d(un1, un2):
                d = bin_dist[uu]
                in1 = (d in bin_dist[un1]) and (-d in bin_dist[un1])
                in2 = (d in bin_dist[un2]) and (-d in bin_dist[un2])
                if in1 and in2:
                    bins_balance += [uu]

            # print(bin_dist[bins_balance])

            # checking that the bins we are using represent each category equally
            assert(np.mean(bin_dist[bins_balance]<0)==0.5)

            n_samp_eachbin = int(np.ceil(n_boot_samp/len(bins_balance)))

            # loop over correct/incorrect trials
            for ci, inds in enumerate([inds1, inds2]):

                nt = len(inds)

                for bi in range(n_boot_iter):

                    # make a resampling order that represents each bin equally
                    inds_resamp = []
                    for bn in bins_balance:
                        inds_bin = inds[coord_binned[inds]==bn]
                        assert(len(inds_bin)>0)
                        # if bi==0:
                        #     print(len(inds_bin), n_samp_eachbin)
                        inds_resamp.append(np.random.choice(inds_bin, n_samp_eachbin, replace=True))    
                    inds_resamp = np.concatenate(inds_resamp, axis=0)
                    # print(len(inds_resamp))

                    # check that the set we created has half each category
                    assert(np.mean(categ_actual[inds_resamp]==1)==0.5)

                    # double check resample order
                    assert(np.all(np.isin(coord_binned[inds_resamp], bins_balance)))
                    counts = np.array([np.sum(coord_binned[inds_resamp]==bn) for bn in bins_balance])
                    assert(np.all(counts==n_samp_eachbin))

                    # get predictions from each ROI, these trials
                    for ri in range(n_rois):

                        pred = dec_withintask['preds_all'][si][ri][ti].astype(int)

                        # figure out what points these 16 values correspond to
                        coords_pred = grid_pts[pred,:].round(2)

                        # binarize the predictions of 16-way classifier into 2 categories
                        # based on current axis "ii"
                        categ_pred = grid_utils.get_categ(coords_pred, (ii+1))

                        prob = dec_withintask['probs_all'][si][ri][ti]

                        # "confidence" in assignment to category 2 vs 1
                        # group the 16 points into categories w/r/t relevant axis
                        g1 = grid_utils.get_categ(grid_pts, ii+1)==1
                        p_categ1 = np.sum(prob[:,g1], axis=1)
                        g2 = grid_utils.get_categ(grid_pts, ii+1)==2
                        p_categ2 = np.sum(prob[:,g2], axis=1)

                        # signed confidence will be: p(correct) - p(incorrect)
                        signedconf = np.zeros_like(p_categ1)
                        signedconf[categ_actual==1] = p_categ1[categ_actual==1] - p_categ2[categ_actual==1]
                        signedconf[categ_actual==2] = p_categ2[categ_actual==2] - p_categ1[categ_actual==2]

                        d = stats_utils.get_dprime(categ_pred[inds_resamp], categ_actual[inds_resamp])
                        dprime_hardtrials_sepcorrect_boot[si,ri,ti,ci,bi] = d

                        signedconf_hardtrials_sepcorrect_boot[si,ri,ti,ci,bi] = np.mean(signedconf[inds_resamp])
    
    # save this, just because it takes a long time to run
    fn2save = os.path.join(save_folder, 'decode_multiclass_sepcorrect_bootstrap.npy')

    np.save(fn2save, {'signedconf_hardtrials_sepcorrect_boot': signedconf_hardtrials_sepcorrect_boot, \
                      'dprime_hardtrials_sepcorrect_boot': dprime_hardtrials_sepcorrect_boot, \
                     })
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pandas as pd

root = '/usr/local/serenceslab/maggie/shapeDim/'

sys.path.append(os.path.join(root, 'Analysis'))
from code_utils import file_utils


def preproc_main_task(sublist = np.arange(1,8)):
    
    subinits = ['S%02.f'%ss for ss in sublist]
    nSubj = len(sublist)

    # % list the orders in which things will happen across entire experiment
    # % [nParts x nSess]
    # % which boundary is active? this defines which categorization dim they're using
    bound_list_all = np.array([[1,2,3],[2,3,1],[3,1,2],[1,2,3],[2,3,1],[3,1,2]]);
    # % which response mapping? e.g. which finger for which category?
    map_list_all = np.array([[1,2,1], [1,2,2], [1,2,1], [2,1,2], [2,1,1], [2,1,2]]);
    # the above are just used to double-check counterbalancing

    # % names of the tasks
    tasklist = ['Linear (1)','Linear (2)','Checker'];
    nTasks = len(tasklist);

    # % info about main categorization tasks
    nRunsPerPart = 2;
    nPartsPerTask = 2;  
    # 2 mappings per task
    nPartsTotal  = nTasks*nPartsPerTask;
    
    possible_responses = (1,2);

    start = 0;    
    stop = 5; 
    center = (stop-start)/2+start;

    # %% quality control settings...
    # % if they fail to respond to 25 percent, probably they were confused and/or
    # % distracted and/or unmotivated.
    # keep these runs in for now, but print a warning
    timeout_pct_cutoff = 0.25;


    for si, ss in enumerate(sublist):

        bdat = pd.DataFrame(columns=['sub','sess','part', 'run_in_part','run_overall',\
                                     'task','map','run_difficulty', \
                                     'trial_in_run','trial_overall', \
                                     'rt','correct_resp','resp','timeout', 
                                     'category_unmapped', 'dist_from_bound_signed', 'resp_unmapped', \
                                     'is_main_grid', 'ptx','pty','quadrant'])

        # figure out which counter-balance condition this subject was in
        cb_ind = np.mod(ss,3);
        
        # get expected counter-balancing order for this subject
        bound_list = bound_list_all
        map_list = map_list_all
        if cb_ind==2:
            # second session goes first
            bound_list = bound_list[:,[1,2,0]];
            map_list = map_list[:,[1,2,0]];
        elif cb_ind==0:
            # third session goes first
            bound_list = bound_list[:,[2,0,1]];
            map_list = map_list[:,[2,0,1]];

        subinit = subinits[si]
        print('\n%s: cb=%d\n'%(subinit, cb_ind))

        # % these are files that i made manually, listing which runs we need to
        # % load for the main task. For some subjects, we're loading a run from
        # % session 3 in place of one from session 1, etc. so this is addressed
        # % in this file...
        main_run_sequence = pd.read_csv(os.path.join(root, 'Samples','run_sequences','%s_main_run_sequence.csv'%subinit))

        # where is the behav data?
        behav_data_folder = os.path.join(root, 'DataBehavior', subinit)

        sess2do = [0,1,2]

        rc = -1; tc=-1 # counters over runs, trials

        for ses in sess2do:

            # looping over the 6 parts of this session
            for xx in np.arange(nPartsTotal):

                # figure out where data for this "part" lives
                ind = np.where((main_run_sequence['sess']==(ses+1)) & \
                               (main_run_sequence['part']==(xx+1)))[0]
                if len(ind)==0:
                    print('Subject %s is missing sess %d part %d\n'%(p['SubNumStr'],ses+1,xx+1))
                    continue
                else:
                    ind = ind[0]
                    sess_date = np.array(main_run_sequence['date'])[ind]
                    actual_sess_num = np.array(main_run_sequence['actual_sess'])[ind]
                    filename = os.path.join(behav_data_folder, 'Session%d'%actual_sess_num, \
                                            '%s_MainTaskMRI_scannerversion_sess%d_part%d_%s.mat'%\
                                            (subinit,(ses+1),xx+1,sess_date));
                    if not os.path.exists(filename):
                        print('Subject %s is missing sess %d part %d\n'%(p['SubNumStr'],ses+1,xx+1))
                        continue
                
                print('loading from %s'%filename)
                
                TheData = file_utils.load_mat_behav_data(filename)
                p = TheData[0]['p']

                # properties of this part
                curr_task = bound_list[xx,ses]
                curr_map = map_list[xx,ses]
                
                # check and make sure the parameters in the file are correct
                assert(p['which_bound']==curr_task);
                assert(p['which_mapping']==curr_map)
                
                run_difficulty = p['RunDifficulty']

                nTrials = len(TheData[0]['data']['Response'])

                for rr in range(nRunsPerPart):

                    rc+=1

                    # check if run is missing (or if we're intentionally skipping, for one run)
                    if (rr>=len(TheData)) or (subinit=='S06' and ses==2 and xx==4 and rr==1):
                        print('Subject %s is missing run %d for sess %d part %d\n'%(p['SubNumStr'],rr+1,ses+1,xx+1))
                        continue

                    p = TheData[rr]['p'];
                    data = TheData[rr]['data'];
                    t = TheData[rr]['t'];
                    
                    response = np.array(data['Response']).astype(float)
                    resptime = np.array(t['RespTimeFromOnset'])
                    
                    # # print(subinit, ses, xx, rr)
                    # if (subinit=='S07') & (ses==2) & (xx==4) & (rr==1):
                    #     # this is a run where the subject had messed up response mapping and this was noted.
                    #     # flip the responses back now
                    #     print('fixing resp mapping')
                    #     response_new = np.full(fill_value=np.nan, shape=np.shape(response))
                    #     response_new[response==2] = 1
                    #     response_new[response==0] = 2
                    #     response = response_new
                    # else:
                    response[np.isnan(resptime)] = np.nan
                    
                    # find trials where subject failed to respond
                    timeout = np.isnan(response) | (response==0) | ~np.isin(response,possible_responses);
                    
                    response[timeout] = np.nan
                    
                    # make sure we didn't mess anything up here
                    assert(not np.any(np.isnan(response[~timeout])))
                    assert(np.all(np.isnan(response[timeout])))
  
                    if np.sum(timeout)>nTrials*timeout_pct_cutoff:
                        # print a warning if there were lots of timeouts
                        print('WARNING: %s sess %d part %d run %d had %d timeout trials'%(subinit,ses+1,xx+1,rr+1,np.sum(timeout)));
                        # print('skipping this run')
                        # continue
                        
                    for tr in range(nTrials):

                        tc+=1
                        rt = resptime[tr]
                        correct_resp = p['category'][tr]
                        resp = response[tr]

                        is_main_grid = p['is_main_grid'][tr]
                        
                        ptx = p['points'][tr][0]
                        pty = p['points'][tr][1]
                        center = 2.5 # define quadrant for each grid pt
                        if (ptx>center) & (pty>center):
                            quadrant = 1;
                        elif (ptx<center) & (pty>center):
                            quadrant = 2
                        elif (ptx<center) & (pty<center):
                            quadrant = 3
                        elif (ptx>center) & (pty<center):
                            quadrant = 4
                        assert(p['quadrant'][tr]==quadrant)
                            
                        # define which category was which - for mapping 2, this
                        # would have been reversed relative to mapping 1, so make
                        # them the same again.
                        if (curr_map)==2:
                            category_unmapped = 3-p['category'][tr]
                        else:
                            category_unmapped = p['category'][tr]

                        if curr_task<3:
                            sign = np.int(np.array(p['points'])[tr,curr_task-1]>center);
                            if sign==0:
                                sign = -1;
                            dist_from_bound = np.round(p['dist_from_bound'][tr],1);
                        else:
                            sign = np.nan; dist_from_bound = np.nan

                        # resp_cat is which category their response
                        # indicated, which is swapped on different
                        # mappings. Make them the same again.
                        if ~timeout[tr]:
                            if (curr_map)==1:
                                resp_unmapped = resp;
                            else:
                                resp_unmapped = 3-resp;
                        else:
                            resp_unmapped = np.nan

                        bdat = bdat.append(pd.DataFrame({'sub': ss, \
                                                         'sess': ses+1, \
                                                         'part':xx+1, \
                                                         'run_in_part': rr+1, \
                                                        'run_overall': rc+1, \
                                                        'task': curr_task, \
                                                        'map': curr_map, \
                                                        'run_difficulty': run_difficulty, \
                                                        'trial_in_run': tr+1, \
                                                        'trial_overall': np.int(tc+1), \
                                                        'rt': rt, \
                                                        'correct_resp': correct_resp, \
                                                        'resp': resp, \
                                                        'timeout': timeout[tr], \
                                                        'category_unmapped': category_unmapped,
                                                        'dist_from_bound_signed': dist_from_bound*sign, \
                                                        'resp_unmapped': resp_unmapped, \
                                                        'is_main_grid': is_main_grid, \
                                                        'ptx': ptx, \
                                                        'pty': pty, \
                                                        'quadrant': quadrant}, \
                                                        index=[tc]))
        fn2save = os.path.join(behav_data_folder, '%s_maintask_preproc_all.csv'%(subinit))
        print('writing to %s'%fn2save)
        bdat.to_csv(fn2save)
        
        
        
        
def preproc_repeat_task(sublist = np.arange(1,8)):

    # sublist = np.arange(1,8)
    subinits = ['S%02.f'%ss for ss in sublist]
    nSubj = len(sublist)

    # %% quality control settings
    # % if they fail to respond to 25 percent, probably they were confused and/or
    # % distracted and/or unmotivated.
    # i'm keeping these runs in for now, but print a warning...
    timeout_pct_cutoff = 0.25;

    possible_responses = [1,2]

    for si, ss in enumerate(sublist):

        bdat = pd.DataFrame(columns=['sub','sess',\
                                     'run_in_sess', 'run_overall',\
                                     'map','run_difficulty', \
                                     'trial_in_run','trial_overall', \
                                     'rt','correct_resp','resp','timeout', 
                                     'is_repeat', \
                                     'dist_from_previous',\
                                     'is_main_grid', \
                                     'ptx','pty', 'quadrant'])
        subinit = subinits[si]
        print('\n%s\n'%(subinit))

        # where is the behav data?
        behav_data_folder = os.path.join(root, 'DataBehavior', subinit)

        sess2do = [0,1,2]

        rc = -1; tc=-1 # counters over runs, trials

        for ses in sess2do:

            # find the behavioral file for this session
            dat_files = os.listdir(os.path.join(behav_data_folder, 'Session%d'%(ses+1)))
            dat_files = [f for f in dat_files if 'OneBackTask' in f]
            assert(len(dat_files)==1)
            filename = os.path.join(behav_data_folder, 'Session%d'%(ses+1), dat_files[0])

            print('loading from %s'%filename)

            TheData = file_utils.load_mat_behav_data(filename)
            n_runs = len(TheData)
            print('found %d runs'%n_runs)
            
            if (subinit=='S06') and (ses==2):
                print('skipping last run in this session')
                # corresponding nifti file is missing for this run,
                # so removing the behavioral file as well.
                n_runs = 4;

            for rr in range(n_runs):

                rc+=1

                p = TheData[rr]['p'];
                data = TheData[rr]['data'];
                t = TheData[rr]['t'];

                curr_map = p['RespMap']
                run_difficulty = p['RunDifficulty']

                nTrials = len(data['Response'])

                response = np.array(data['Response']).astype(float)
                resptime = np.array(t['RespTimeFromOnset'])

                response[np.isnan(resptime)] = np.nan

                # find trials where subject failed to respond
                timeout = np.isnan(response) | (response==0) | ~np.isin(response,possible_responses);

                response[timeout] = np.nan

                # make sure we didn't mess anything up here
                assert(not np.any(np.isnan(response[~timeout])))
                assert(np.all(np.isnan(response[timeout])))

                if np.sum(timeout)>nTrials*timeout_pct_cutoff:
                    # print a warning if there were lots of timeouts
                    print('WARNING: %s sess %d run %d had %d timeout trials'%(subinit,ses+1,rr+1,np.sum(timeout)));
                    # print('skipping this run')
                    # continue

                for tr in range(nTrials):

                    tc+=1
                    rt = resptime[tr]
                    correct_resp = p['correct_resp'][tr]
                    is_repeat = p['is_repeat'][tr]
                    resp = response[tr]

                    is_main_grid = p['is_main_grid'][tr]

                    ptx = p['points'][tr][0]
                    pty = p['points'][tr][1]
                    center = 2.5 # define quadrant for each grid pt
                    if (ptx>center) & (pty>center):
                        quadrant = 1;
                    elif (ptx<center) & (pty>center):
                        quadrant = 2
                    elif (ptx<center) & (pty<center):
                        quadrant = 3
                    elif (ptx>center) & (pty<center):
                        quadrant = 4
                    
                    dist_from_previous = p['dist_from_previous'][tr]

                    bdat = bdat.append(pd.DataFrame({'sub': ss, \
                                                    'sess': ses+1, \
                                                    'run_in_sess': rr+1, \
                                                    'run_overall': rc+1, \
                                                    'map': curr_map, \
                                                    'run_difficulty': run_difficulty, \
                                                    'trial_in_run': tr+1, \
                                                    'trial_overall': np.int(tc+1), \
                                                    'rt': rt, \
                                                    'correct_resp': correct_resp, \
                                                    'resp': resp, \
                                                    'timeout': timeout[tr], \
                                                    'is_repeat': is_repeat, \
                                                    'dist_from_previous': dist_from_previous, \
                                                    'is_main_grid': is_main_grid, \
                                                    'ptx': ptx, \
                                                    'pty': pty, \
                                                    'quadrant': quadrant}, \
                                                    index=[tc]))

        fn2save = os.path.join(behav_data_folder, '%s_reptask_preproc_all.csv'%(subinit))
        print('writing to %s'%fn2save)
        bdat.to_csv(fn2save)
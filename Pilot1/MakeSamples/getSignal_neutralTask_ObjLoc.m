% extract signal in each voxel averaging over several TRs following trial
% events, save values for input to further analysis

clear
close all;
root = '/usr/local/serenceslab/maggie/shapeDim/Pilot1/Samples/';

sublist = {'01'};

for ss=1:length(sublist)

    mainSig = struct('trialLabs',[]);

    my_subject = sublist{ss};

    load([root, 'SampleFile_ObjLoc_S', my_subject,'.mat'],'samplesMain','main','ROIs','all_vox_concat');

    fn2save = [root 'NeutralTaskSignalByTrial_ObjLoc_S' my_subject, '.mat'];
   
    %% vis areas to use

    my_areas = {'fusiform','lateraloccipital','inferiorparietal','superiorparietal'};
   
    trDur = .8; % actual duration of the TR, for plotting, in seconds...
    nTRs = 387 - 16;   

    %% how many TRs to average?
    % 5 to 7 TRs after event onset, 
    avgTRs_stim = [5,7];
    %%
    for vv = 1:length(my_areas) % for all visual areas I want to look at
        %% pull out the data from each ROI
        % row is the hemisphere
        [row_inds,col_inds] = find(reshape(contains({ROIs.name}, my_areas(vv)),2,[]));
        fullDat=[]; 
        
        for ii=1:length(row_inds)
            name = ROIs(row_inds(ii),col_inds(ii)).name;
            if ~isempty(ROIs(row_inds(ii),col_inds(ii)).voxel_inds)
                % jj gives indices into the all_vox_concat array
                [~,jj]=intersect(all_vox_concat, ROIs(row_inds(ii),col_inds(ii)).voxel_inds);
                
                fullDat = [fullDat, samplesMain(:,jj)];               
            end
        end
        
        nVox = size(fullDat,2);

        if nVox==0
            fprintf('no voxels in area %s!\n',my_areas{vv});
            continue
        end

        fprintf('processing area %s, %d voxels\n', my_areas{vv}, nVox);
        
        %% now zscore the data from each run to normalize...

        nRuns = size(fullDat,1)/nTRs; 
        
        scan_nums = zeros(size(fullDat,1),1);
    
        if mod(nRuns,1)~=0
            error('something bad happened here with fullDat run length')
        end
        for ii=1:nRuns
            fullDat(ii*nTRs-nTRs+1:ii*nTRs,:) = zscore(fullDat(ii*nTRs-nTRs+1:ii*nTRs, :));
            scan_nums(ii*nTRs-nTRs+1:ii*nTRs,:) = ii;
        end

       %% label the data

        % find the start of each trial (grating onset) - will be marked by
        % a shift from label 0 to 1
        % first reshape this list of event labels [nTRs x nRuns]
        event_labels_reshaped = reshape(main.EventLabels,nTRs,length(main.EventLabels)/nTRs);
        % put a row of zeros before it and get a diff - this ensures the
        % first event gets detected.
        trial_onset_bool = diff([zeros(1,size(event_labels_reshaped,2)); event_labels_reshaped],[],1)==1;
        trial_onset_bool = trial_onset_bool(:);
        
        % list which TR these events occured at - use this to get trial
        % labels.
        trial_onset_num = find(trial_onset_bool(:));  
      
        %% save out a bunch of descriptors for the trials        
         
        mainSig(vv).StimLabels = main.StimLabels(trial_onset_num,:);
        mainSig(vv).ShapeCoords = main.ShapeCoords(trial_onset_num,:);
       
        mainSig(vv).isRepeat = main.RepeatLabels(trial_onset_num,:);

        mainSig(vv).runLabs = main.RunLabels(trial_onset_num);
        mainSig(vv).trialLabs = main.TrialLabels(trial_onset_num);

        mainSig(vv).scanLabs = scan_nums;

        %% avg the data across each trial
        % We do this to come up with 1 response value for each voxel on each 
        % trial - this will be input to the IEM.
        % Do this on a run-by-run basis so that we can make sure we're not 
        % grabbing timepoints in our avg that extend beyond the end of the run. 
        % Not the fastest routine, but ok...
        nTrials = numel(trial_onset_num);

        assert(nTrials==48*16)
        
        mainDat = nan(nTrials, nVox); 
       
        triCnt = 0; % counter across "good" trials where the entire desired avg window is available.
    
        for rr=1:nRuns

            runInd = rr*nTRs-nTRs+1:rr*nTRs;
            assert(all(find(main.RunLabels==rr)==runInd'))
            curDat = fullDat(runInd,:);      % data from current run

            % get numeric indices for each event in this run
            these_trial_onsets = find(trial_onset_bool(runInd));       
           
            
            for tt=1:numel(these_trial_onsets)
                % check to make sure we don't go past the end of the run.
                if these_trial_onsets(tt)+avgTRs_stim(2)<=nTRs
                    triCnt = triCnt + 1;  % increment counter over good trials

                    % sum the data over desired interval.
                    mainDat(triCnt,:) = mean(curDat(these_trial_onsets(tt)+avgTRs_stim(1):these_trial_onsets(tt)+avgTRs_stim(2),:));
                    
                end
            end
        end

        assert(triCnt==nTrials)
        assert(~sum(isnan(mainDat(:))))

        mainSig(vv).mainDat = mainDat;
       
    end

    fprintf('saving to %s\n',fn2save);
    save(fn2save,'mainSig','my_areas');

end
% Get the timing of all events in all runs
% saves out a file that is used by MakeSampleFile and by RunGLM_WholeBrain

% MMH 9/5/18

clear
close all

% set inputs
FS_sbj = 'AV';
subnum = '01';

% set paths (end all with filesep!)
experiment_path = '/mnt/neurocube/local/serenceslab/maggie/faceDim/pilot4/';
out_path = [experiment_path, 'Samples/'];
beh_path = [experiment_path, 'DataBehavior/'];

%% Timing information

TRdur = .8;
TRditched = 16; % 16 TRs were removed
not_recorded = TRdur * TRditched;

%% Main Task

nRuns = 0; nTRs_main = 319 - TRditched; 
RunLabels = []; EventLabels = []; FaceLabels = []; TrialLabels = []; RTLabels = [];
RepeatLabels = [];
ResponseLabels = [];

trialnumglobal = 0;

% get data file from each session (they are not concatenated here)
loc_file = dir([beh_path, 'S' char(subnum) '/Session*/*FDA_MatchTask*']);
actual_sess = [1,2,3];
for sess = 1:length(loc_file)
    
load([beh_path 'S' char(subnum) filesep 'Session' num2str(actual_sess(sess)) filesep loc_file(sess).name])

if sess==1
    nTrialsEach = TheData(1).p.NumTrials;
end

for run = 1:length(TheData)

    nRuns = nRuns + 1;

    % Find stimulus onset time
    Trial_onset = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime  -not_recorded;
    Trial_offset = TheData(run).t.stim_flips(:,6) - TheData(run).t.StartTime - not_recorded;
    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;
    
%     theseglobalinds =[];
    for n = 1:TheData(run).p.NumTrials
        
        trialnumglobal = trialnumglobal+1;
%         theseglobalinds = [ theseglobalinds trialnumglobal];
%         
%         if ismember(n,TheData(run).p.trials2switch)
%             
%             mbcount = find(TheData(run).p.trials2switch==n);
            
            % list times for all events on this trial - including the
            % pre-trial interblock interval if there is one
            event_times_this_trial = [Trial_onset(n),Trial_offset(n)];
            
            Event_onset = [Event_onset, event_times_this_trial];
            % n.1 is miniblock cue, 1 is grating, 2 is dial
            Event_type = [Event_type, 1, 0]; 
            % we will mark the miniblock cue interval with the trial
            % number+0.1. This will discriminate it from the actual trial
            % onset, which we need down below.
            % also mark the dial onset time with n-0.1 - this will allow us
            % to make sure that every event gets marked, even if they are
            % sub-TR length!
            Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,2)];  
            Event_trial_local = [Event_trial_local, repmat(n,1,2)];

    end
    
    % Find onset times for all my TRs
    TR_onset = 0:.8:nTRs_main*0.8;
    TR_onset = TR_onset(1:end-1);

    % list everything by TRs
    triallabslocal_byTR = zeros(length(TR_onset),1);
    eventlabs_byTR = zeros(length(TR_onset),1);
    triallabsglobal_byTR = zeros(length(TR_onset),1);
    
    for i = 1:length(TR_onset)
        
        middleofTR = TR_onset(i)+(TRdur/2);
        myind = find(middleofTR>Event_onset, 1, 'last');

        % on this TR - what type of of event is happening, and what global
        % trial number and run number is it a part of?
        eventlabs_byTR(i) = Event_type(myind);
        triallabslocal_byTR(i) = Event_trial_local(myind);
        triallabsglobal_byTR(i) = Event_trial_global(myind);
        
    end
    
    % make sure that there is a 1 marking the start of every image. This is necessary
    % because sometimes the events might not land in the right place in a
    % TR, because they may be <0.8 s duration.
    for n = 1:nTrialsEach           
        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
    end
   
      
    % now round all the trial numbers back to being integers (need this for
    % the code below where we use them as indices)
    triallabslocal_byTR = round(triallabslocal_byTR,0);
    triallabsglobal_byTR = round(triallabsglobal_byTR,0);
    
    numZeros = sum(triallabslocal_byTR==0);
    
    RunLabels = [RunLabels;repmat(nRuns,length(TR_onset),1)];
    TrialLabels = [TrialLabels;triallabsglobal_byTR];
    EventLabels = [EventLabels; eventlabs_byTR];

    RepeatLabels = [RepeatLabels;  nan(numZeros,1); TheData(run).p.isRepeat(triallabslocal_byTR(numZeros+1:end))];
    
    FaceLabels = [FaceLabels;  nan(numZeros,1); TheData(run).p.faceList(triallabslocal_byTR(numZeros+1:end))];
    RTLabels = [RTLabels; nan(numZeros,1); TheData(run).t.RespTimeFromOnset(triallabslocal_byTR(numZeros+1:end))];
    ResponseLabels = [ResponseLabels; nan(numZeros,1); TheData(run).data.Response(triallabslocal_byTR(numZeros+1:end))];
    
end %run loop

end

main.RunLabels = RunLabels;
main.TrialLabels = TrialLabels;
main.EventLabels = EventLabels;
main.FaceLabels = FaceLabels;
main.RTLabels = RTLabels;
main.RepeatLabels = RepeatLabels;
main.ResponseLabels = ResponseLabels;

% convert these into multidimensional feature vectors
main.FaceNumOrig = TheData(1).p.faceNums2Use(FaceLabels,:);
main.PolCoords = TheData(1).p.pol_coords(FaceLabels,:);
main.RaceGenderCoords = TheData(1).p.cart_coords(FaceLabels,:);

if numel(main.EventLabels)~=nTRs_main*numel(unique(main.RunLabels))
    error('wrong number of total TRs!')
end

%% ATTention tasks

% nRuns = 0;
nTRs_att = 402 - TRditched; 
RunLabels = []; EventLabels = []; FaceLabels = []; TrialLabels = []; RTLabels = [];
TaskLabels = [];
ResponseLabels = [];
RespMapLabels = [];
CorrectRespLabels = [];
SessLabels = [];

trialnumglobal = 0;

% get data file from each session (they are not concatenated here)
loc_file = dir([beh_path, 'S' char(subnum) '/Session*/*FDA_AttnTask*']);
actual_sess = [1,2,3];
for sess = 1:length(loc_file)
    
load([beh_path 'S' char(subnum) filesep 'Session' num2str(actual_sess(sess)) filesep loc_file(sess).name])

if sess==1
    nTrialsEach = TheData(1).p.NumTrials;
    nRunPairsFirstSess = length(TheData)/2;
end

for run = 1:length(TheData)
    
    if sess==1
        runNumber = TheData(run).p.RunPairNum;
    elseif sess==2
        runNumber = nRunPairsFirstSess+TheData(run).p.RunPairNum;
    elseif sess==3
        runNumber = nRunPairsFirstSess*2+TheData(run).p.RunPairNum;
    end

%     nRuns = nRuns + 1;

    % Find stimulus onset time
    Trial_onset = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime  -not_recorded;
    Trial_offset = TheData(run).t.stim_flips(:,6) - TheData(run).t.StartTime - not_recorded;
    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;
    
%     theseglobalinds =[];
    for n = 1:TheData(run).p.NumTrials
        
        trialnumglobal = trialnumglobal+1;
%         theseglobalinds = [ theseglobalinds trialnumglobal];
%         
%         if ismember(n,TheData(run).p.trials2switch)
%             
%             mbcount = find(TheData(run).p.trials2switch==n);
            
            % list times for all events on this trial - including the
            % pre-trial interblock interval if there is one
            event_times_this_trial = [Trial_onset(n),Trial_offset(n)];
            
            Event_onset = [Event_onset, event_times_this_trial];
            % n.1 is miniblock cue, 1 is grating, 2 is dial
            Event_type = [Event_type, 1, 0]; 
           
            Event_trial_global = [Event_trial_global, repmat(trialnumglobal,1,2)];  
            Event_trial_local = [Event_trial_local, repmat(n,1,2)];

    end
    
    % Find onset times for all my TRs
    TR_onset = 0:.8:nTRs_att*0.8;
    TR_onset = TR_onset(1:end-1);

    % list everything by TRs
    triallabslocal_byTR = zeros(length(TR_onset),1);
    eventlabs_byTR = zeros(length(TR_onset),1);
    triallabsglobal_byTR = zeros(length(TR_onset),1);
    
    for i = 1:length(TR_onset)
        
        middleofTR = TR_onset(i)+(TRdur/2);
        myind = find(middleofTR>Event_onset, 1, 'last');

        % on this TR - what type of of event is happening, and what global
        % trial number and run number is it a part of?
        eventlabs_byTR(i) = Event_type(myind);
        triallabslocal_byTR(i) = Event_trial_local(myind);
        triallabsglobal_byTR(i) = Event_trial_global(myind);
        
    end
    
    % make sure that there is a 1 marking the start of every image. This is necessary
    % because sometimes the events might not land in the right place in a
    % TR, because they may be <0.8 s duration.
    for n = 1:nTrialsEach           
        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
    end
   
      
    % now round all the trial numbers back to being integers (need this for
    % the code below where we use them as indices)
    triallabslocal_byTR = round(triallabslocal_byTR,0);
    triallabsglobal_byTR = round(triallabsglobal_byTR,0);
    
    numZeros = sum(triallabslocal_byTR==0);
    
    RunLabels = [RunLabels;repmat(runNumber,length(TR_onset),1)];
    SessLabels = [SessLabels;repmat(sess, length(TR_onset), 1)];
    RespMapLabels = [RespMapLabels; repmat(TheData(run).p.respMap,length(TR_onset),1)];
    TaskLabels = [TaskLabels; repmat(TheData(run).p.TaskNum, length(TR_onset),1)];
    
    if TheData(run).p.TaskNum==1 && (~contains(TheData(run).p.InstrText, 'Race') || contains(TheData(run).p.InstrText, 'Gender'))
        error('task number is labeled wrong')
    end
    if TheData(run).p.TaskNum==2 && (~contains(TheData(run).p.InstrText, 'Gender') || contains(TheData(run).p.InstrText, 'Race'))
        error('task number is labeled wrong')
    end
    
    TrialLabels = [TrialLabels;triallabsglobal_byTR];
    EventLabels = [EventLabels; eventlabs_byTR];

    FaceLabels = [FaceLabels;  nan(numZeros,1); TheData(run).p.faceList(triallabslocal_byTR(numZeros+1:end))];
    RTLabels = [RTLabels; nan(numZeros,1); TheData(run).t.RespTimeFromOnset(triallabslocal_byTR(numZeros+1:end))];
    ResponseLabels = [ResponseLabels; nan(numZeros,1); TheData(run).data.Response(triallabslocal_byTR(numZeros+1:end))];
    
    CorrectRespLabels = [CorrectRespLabels; nan(numZeros,1); TheData(run).p.correctResp(triallabslocal_byTR(numZeros+1:end))'];
    
end %run loop

end


assert(all(unique(SessLabels)'==1:3))
assert(all(unique(RunLabels)'==1:18))

att.SessLabels = SessLabels;
att.RunLabels = RunLabels;

att.TrialLabels = TrialLabels;
att.EventLabels = EventLabels;
att.FaceLabels = FaceLabels;
att.RTLabels = RTLabels;
att.ResponseLabels = ResponseLabels;
att.TaskLabels = TaskLabels;
att.RespMapLabels = RespMapLabels;
att.CorrectRespLabels =CorrectRespLabels;

% convert these into multidimensional feature vectors
% in the trialInfo.featCombs matrix - [1=negative, 2=positive, 3=zero]
att.FaceNumOrig = trialInfo.faceNums2Use(FaceLabels,:);
att.PolCoords = trialInfo.pol_coords(FaceLabels,:);
att.RaceGenderCoords = trialInfo.race_gender_coords(FaceLabels,:);
att.CategoryLabs = trialInfo.category_labs(FaceLabels,:);

% check and make sure categories are all correct
expLabs = zeros(size(att.CategoryLabs,1),1);
expLabs(att.TaskLabels==1) = att.CategoryLabs(att.TaskLabels==1,1);
expLabs(att.TaskLabels==2) = att.CategoryLabs(att.TaskLabels==2,2);
expLabs(expLabs~=0 & att.RespMapLabels==2) = 3-expLabs(expLabs~=0 & att.RespMapLabels==2);

same = expLabs==att.CorrectRespLabels;
if any(~same(expLabs~=0))
    error('error assigning category labels')
end

if numel(att.EventLabels)~=nTRs_att*trialnumglobal/nTrialsEach
    error('wrong number of total TRs!')
end
%% Save timing file
filename = ['TimingFile_S', subnum];
fprintf('saving file to %s\n',filename);

if ~exist(out_path, 'dir'), mkdir(out_path); end
save([out_path, filename], 'main', 'att','-v7.3');


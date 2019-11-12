% Get the timing of all events in all runs
% saves out a file that is used by MakeSampleFile and by RunGLM_WholeBrain

% MMH 9/5/18

clear
close all

% set inputs
FS_sbj = 'CG';
subnum = '01';

% set paths (end all with filesep!)
experiment_path = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Pilot1/';
out_path = [experiment_path, 'Samples/'];
beh_path = [experiment_path, 'DataBehavior/'];

%% Timing information

TRdur = .8;
TRditched = 16; % 16 TRs were removed
not_recorded = TRdur * TRditched;

%% Main Task

nRuns = 0; nTRs_main = 387 - TRditched; 
RunLabels = []; EventLabels = []; StimLabels = []; TrialLabels = []; RTLabels = [];
RepeatLabels = [];
ResponseLabels = [];

trialnumglobal = 0;

% get data file from each session (they are not concatenated here)
loc_file = dir([beh_path, 'S' char(subnum) '/Session*/*OneBackPilotTask*']);
actual_sess = [1];
for sess = 1:length(loc_file)
    
load([beh_path 'S' char(subnum) filesep 'Session' num2str(actual_sess(sess)) filesep loc_file(sess).name])

if sess==1
    nTrialsEach = TheData(1).p.nTrials;
end

for run = 1:length(TheData)

    nRuns = nRuns + 1;

    % Find stimulus onset time
    Trial_onset = TheData(run).t.stim_flips(:,1) - TheData(run).t.StartTime  -not_recorded;
    Trial_offset = TheData(run).t.stim_flips(:,2) - TheData(run).t.StartTime - not_recorded;
    Event_onset = 0;  Event_type = 0; Event_trial_global = 0; Event_trial_local = 0;
 
    for n = 1:TheData(run).p.nTrials
        
        trialnumglobal = trialnumglobal+1;

        % list times for all events on this trial 
        event_times_this_trial = [Trial_onset(n),Trial_offset(n)];

        Event_onset = [Event_onset, event_times_this_trial];
        % 1 is stim on, 0 is stim off
        Event_type = [Event_type, 1, 0]; 
        % also appending to a list, that is nEvents long and describes which trial each event is from 
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
    % TR
    for n = 1:nTrialsEach           
        eventlabs_byTR(find(triallabslocal_byTR==n,1))= 1;
    end
   
      
    % now round all the trial numbers back to being integers (need this for
    % the code below where we use them as indices)
    triallabslocal_byTR = round(triallabslocal_byTR,0);
    triallabsglobal_byTR = round(triallabsglobal_byTR,0);
    
    % how many TRs elapse before the first event? use this to keep all my
    % arrays the same length
    numZeros = sum(triallabslocal_byTR==0);
    
    % append to these big lists, which are nTRs long by the end
    RunLabels = [RunLabels;repmat(nRuns,length(TR_onset),1)];
    TrialLabels = [TrialLabels;triallabsglobal_byTR];    
    EventLabels = [EventLabels; eventlabs_byTR];
    
    % here i'm taking out info from my data structure, using my trial
    % indices created above to get from trial space to TR space
    RepeatLabels = [RepeatLabels;  nan(numZeros,1); TheData(run).p.isrepeat(triallabslocal_byTR(numZeros+1:end))];
    StimLabels = [StimLabels;  nan(numZeros,1); TheData(run).p.imlist(triallabslocal_byTR(numZeros+1:end))];
    RTLabels = [RTLabels; nan(numZeros,1); TheData(run).t.RespTimeFromOnset(triallabslocal_byTR(numZeros+1:end))];
    ResponseLabels = [ResponseLabels; nan(numZeros,1); TheData(run).data.Response(triallabslocal_byTR(numZeros+1:end))];
    
end %run loop

end

% now i have a bunch of lists nTRsTotal long, and they will all get put
% into my main data structure for saving. 
main.RunLabels = RunLabels;
main.TrialLabels = TrialLabels;
main.EventLabels = EventLabels;
main.StimLabels = StimLabels;
main.RTLabels = RTLabels;
main.RepeatLabels = RepeatLabels;
main.ResponseLabels = ResponseLabels;

% convert these into multidimensional feature vectors
main.ShapeCoords = TheData(1).p.points(StimLabels,:);

if numel(main.EventLabels)~=nTRs_main*numel(unique(main.RunLabels))
    error('wrong number of total TRs!')
end

%% Save timing file
filename = ['TimingFile_S', subnum];
fprintf('saving file to %s\n',filename);

if ~exist(out_path, 'dir'), mkdir(out_path); end
save([out_path, filename], 'main','-v7.3');


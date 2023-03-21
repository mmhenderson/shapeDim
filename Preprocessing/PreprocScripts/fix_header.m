{\rtf1\ansi\ansicpg1252\cocoartf2513
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 Monaco;}
{\colortbl;\red255\green255\blue255;\red22\green21\blue22;\red22\green21\blue22;}
{\*\expandedcolortbl;;\cssrgb\c11373\c10980\c11373;\cssrgb\c11373\c10980\c11373\c3922;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl360\partightenfactor0

\f0\fs24 \cf2 \cb3 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 clear\
close all\
\
% find my root directory - up a few dirs from where i am now\
mypath = pwd;\
filesepinds = find(mypath==filesep);\
nDirsUp = 2;\
exp_path = mypath(1:filesepinds(end-nDirsUp+1));\
inID = 'topup.'; % will use this string to look for raw functional data (nifti)\
\
sub_list_big = \{'S01'\};\
sub_doreti_list_big = \{'CI'\};\
all_sessions_big = \{[1 2 3]\};   % list all sessions that exist for this sub (get merged to runs.list together)\
preproc_sessions_big = \{[1 2 3]\}; % list the ones you want to do right now (those not done yet)\
\
% pathToCopyFrom = fullfile(exp_path,'DataPreprocTest','S01','01_REG_MC_001.nii.gz');\
% pathToCopyTo1 = fullfile(exp_path,'DataPreprocTest','S01','02_REG_MC_001_copy.nii.gz');\
% pathToCopyTo2 = fullfile(exp_path,'DataPreprocTest','S01','03_REG_MC_001_copy.nii.gz');\
\
%% \
\
for xx = 1:length(sub_list_big)\
    \
    sub = sub_list_big\{xx\};\
    sub_doreti = sub_doreti_list_big\{xx\};\
    all_sessions = all_sessions_big\{xx\};\
    preproc_sessions = preproc_sessions_big\{xx\};\
    \
    subpath = [exp_path, 'DataPreproc/', sub, '/']; % path to this subject's preprocessed data\
    \
    raw_path = cell(length(all_sessions),1); % Will be used to get runslist\
    invols = cell(length(all_sessions),1);\
    runnrs = cell(length(all_sessions),1);\
    \
    pathToCopyHeaderFrom = [subpath, '01_REG_MC_001.nii.gz'];\
    \
    %% iterate through sessions! \
    \
    for session = all_sessions\
        raw_path\{session\} = [exp_path, 'DataRaw/', sub, '/Session'  num2str(session), '/Niftis/']; %path to raw data\
        \
        % Read in the runs.list for this session\
        runnrs_tmp = [];\
        fid = fopen([char(raw_path(session)), 'runs.list']);\
        line = fgetl(fid);\
        while ischar(line) && ~isempty(line)\
            runnrs_tmp = [runnrs_tmp; line(1:3)];\
            line = fgetl(fid);\
        end\
        fclose(fid);\
        runnrs\{session\} = runnrs_tmp;\
        \
        % determine the runs to go into this part\
        clear invols % invols was previously defined as the raw data\
        invols = dir([subpath, sprintf('%02d',session), '_REG_*nii.gz']); % here we will instead use the transformed data\
        % MMH adding this 5/13/19 - check to make sure we're only looking\
        % at runs which are NOT motion corrected yet. Otherwise if we\
        % have a crash in the middle of this script and try to restart it,\
        % we'd get an error.\
        invols = invols(~contains(\{invols.name\},'MC'));\
        if size(runnrs\{session\},1)~=length(invols) % make sure I have the right number of invols because I'm double checking everything like stupid\
            fprintf(['  WARNING: The number of runs in your runs.list does not match the number of transformed nifti''s for session ', num2str(session), '...\\n\\n  RETURNING\\n\\n']);\
            return;\
        end\
        \
\
        for func_idx = 1:length(invols)\
            pathToCopyHeaderTo = ...\
                [subpath, sprintf('%02d',session), '_REG_MC_', runnrs\{session\}(func_idx,:), '.nii.gz']; %open the .par file where the translations and rotations are stored\
                \
        unix(['fslcpgeom ' pathToCopyHeaderFrom ' ' pathToCopyHeaderTo]);\
        fprintf('Copied new header to sess %d run %d\\n',session,func_idx);\
        end\
    end\
    \
end\
}
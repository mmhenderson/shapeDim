function ShapeDim_unwarp_forSiemens(inpath, domoco)
% recon does unwarping of reconstructed multiband data.
% It looks for Niftis in inpath/Niftis and for DICOM in inpath/Dicoms (you
% can change these defaults if you want by changing the strings that are
% entered for nii_dir and dcm_dir line 39 and 40, respectively). The 
% function writes unwarped .nii.gz files to inpath/Niftis. The moco flag is
% off by default, however, if you want the CFMRI routine to do moco (using
% AFNI) for you before the unwarping, you can opt into this by setting this
% flag to 1. Finally, this function also writes a log files to your 
% outpath, specifying the files that went into each unwarping call.
%
% inpath =
% '/mnt/neurocube/local/serenceslab/Rosanne/WM_Distract/DataRaw/S04/Session1'
%

dbstop if error

%% catch errors
% If there is no outpath, write to the inpath
if nargin < 2, domoco = false; end
%%
% If there are no file seperators ('/') at the end of your in or oupaths, put one there:
if ~strcmp(inpath(end), filesep)
    inpath = [inpath filesep];
end

%% Set defaults
my_path = pwd; % I'm here
addpath(fullfile(pwd, 'jsonlab-1.5'))
% topup scripts path. Notice how I have my own copy of all the scripts...
% this is just a precaution. A tmp folder will be made wherever the script
% is ran from, and in this function I cd to the outpath, so a bunch of junk 
% outout should go there. But just in the unlikely case stuff does get 
% written to the path where the topup scrips live, two people doing recon/
% unwarping at the same time, might interfere in unexpected ways. 
topup_dir = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/Preprocessing/PreprocScripts/unwarp_forSiemens/'; 
addpath(genpath(topup_dir));    % topup scripts path
nii_dir = 'Niftis';             % my reconstructed will be found under inpath/nii_dir
dcm_dir = 'Dicoms';             % my dicoms will be found under inpath/dcm_dir
numcores = 8;
if domoco, outstr = '_topup_moco'; mocostr = ' -domoco'; else outstr = '_topup'; mocostr = ''; end

%% Read reconstructed Nifti folder names (functional data)

% Let's get a list of our functional niftis -- 
niis = dir([inpath, nii_dir, '/*serences.nii.gz']);

num_niis = size(niis,1);

%% Read the top-up Nifti directories
% Make sure the FIRST file is the SAME DIRECTION as the main functional runs (posterior to anterior) 
% That is, top_updirs(1) = PA, top_updirs(2) = AP

% KA -- let's make this functional with more than 1 set of top-ups!!! 
% topup_dirs(1) = dir([inpath, nii_dir, '/*DISTORTION_PA.nii.gz']);
% topup_dirs(2) = dir([inpath, nii_dir, '/*DISTORTION_AP.nii.gz']);

% Separate structures for forward and reverse top-up dirs... Make this okay
% if we end up having more than one set!! 
topup_dirs_forward = dir([inpath, nii_dir, '/*DISTORTION_PA.nii.gz']);
topup_dirs_reverse = dir([inpath, nii_dir, '/*DISTORTION_AP.nii.gz']);

% if there is more than one top-up, we'll want to use the one that is
% closest to our functional run
n_topups = length(topup_dirs_forward);
topup_nums = NaN(1,n_topups);
for run = 1:n_topups
    topup_nums(run) = str2num(topup_dirs_forward(run).name(1:4));
end

%% Get the total readout time from the json file
% MAKE SURE TO CHECK THIS IS CORRECT WHEN SETTING UP A NEW PROTOCOL!! 
topup_dir_header =  dir([inpath, nii_dir, '/*DISTORTION_PA.json']);
headerinfo = loadjson([inpath, nii_dir, filesep, topup_dir_header(1).name]);
readout_time = headerinfo.TotalReadoutTime;

%% Check if unwarped files already exist before starting
do_topup = ones(num_niis,1);
for file_idx = 1:num_niis
    if exist([inpath, 'Niftis/', niis(file_idx).name(1:end-7), outstr, '.nii.gz'],'file')
        rsp =  input(['A topupped nifti (',niis(file_idx).name(1:end-7), outstr, '.nii.gz) already exists. Overwrite? (y/n) '], 's');
        if strcmp(rsp, 'n')
            do_topup(file_idx) = 0;
        else
            unix(['rm ', niis(file_idx).name(1:end-7), outstr, '.nii.gz']);
        end
    end
end
    
%% Do unwarping   
topuplogfile = [inpath, 'Niftis/topup.log'];
fid = fopen(topuplogfile,'w+');
fprintf(fid,'Files used for topup\n.nii.gz\t\t\t\tfwd\trvs\n');
for file_idx = 1:num_niis
    if do_topup(file_idx)
        
        % Get the current run number, and decide which topup to use
        run_num = str2num(niis(file_idx).name(1:4));
        [min_diff,min_idx] = min(abs(topup_nums-run_num));
        
        % log the input files to my recon call
        fprintf(fid,['\n', niis(file_idx).name, '\t', topup_dirs_forward(min_idx).name, '\t', topup_dirs_reverse(min_idx).name]);
        % if they don't exist in the outpath yet, move my topup sDirs there
        cd([inpath, 'Niftis']) % <--important, needs all its shit in one place in order to work
        
        % and if a tmp folder is in my outpath, get rid of it
        if exist([inpath, nii_dir, 'tmp'],'dir')
            unix(['rm -rf ', inpath, nii_dir, 'tmp'],'-echo');
        end
        %------------------------------------------------------------------
        % actual unwarping
        %------------------------------------------------------------------
        file_to_unwarp = niis(file_idx).name;
        
        [s,~] = unix([topup_dir, 'run_topup_JW -d1 ',inpath, nii_dir, filesep, topup_dirs_forward(min_idx).name(1:end-7),...
            ' -d2 ', inpath, nii_dir, filesep, topup_dirs_reverse(min_idx).name(1:end-7),...
            ' -i ', inpath, nii_dir, filesep, file_to_unwarp(1:end-7),...
            ' -o ', file_to_unwarp(1:end-7), outstr, mocostr], '-echo');
        %------------------------------------------------------------------
        
        % check for error
        if s > 0 
            fprintf(['\nUnwarping (a.k.a. toptup) failed for ', file_to_unwarp, '.nii.gz\n\n'])
        else % if no error I'm gonna clean up
            stupid_logs = dir(['topup_*.log']);
            if ~isempty(stupid_logs)
                for numlogs = 1:length(stupid_logs)
                    unix(['rm ', stupid_logs(numlogs).name],'-echo');
                end
            end
        end
    end
end
fclose(fid); 
   

%% return to where I started
cd(my_path)
   

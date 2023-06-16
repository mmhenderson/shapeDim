#!/bin/bash
#SBATCH --partition=general
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=1-00:00:00

echo $SLURM_JOBID
echo $SLURM_NODELIST

# activate my conda env here
# need to use full path, this is path to env on SSRDE
source activate /cube/neurocube/local/serenceslab/maggie/conda_envs/shape_dim

# my project path
ROOT=/cube/neurocube/local/serenceslab/maggie/shapeDim/

# put the code directory on your python path
PYTHONPATH=:${ROOT}:Analysis/image_similarity
# /home/AD/mmhender/anaconda3/lib/python3.6/:/home/AD/mmhender/anaconda3/lib/python3.6/site-packages/

# go to folder where script is located
cd ${ROOT}Analysis/image_similarity

debug=0
overwrite=1

# python3 -c 'import image_utils; image_utils.prep_brick('${debug}')'

gist_code_dir=${ROOT}Analysis/image_similarity/gist_matlab/
cd ${gist_code_dir}

matlabpath=/cube/neurocube/local/MATLAB/R2018b/bin/matlab

${matlabpath} -nodisplay -nodesktop -nosplash -r "get_gist($debug, $overwrite); exit"

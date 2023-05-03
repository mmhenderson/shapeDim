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

# change this path
ROOT=/cube/neurocube/local/serenceslab/maggie/shapeDim/

# put the code directory on your python path
PYTHONPATH=:${ROOT}Analysis/

# go to folder where script is located
cd ${ROOT}Analysis/decoding/

debug=0
n_threads=8

python3 -c 'from decoding import decode_binary_withintask; decode_binary_withintask.decode('${debug}', '${n_threads}')'
#!/bin/bash
#SBATCH --partition=general
#SBATCH --gres=gpu:0
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --output=./sbatch_output/output-%A-%x-%u.out 
#SBATCH --time=7-00:00:00

echo $SLURM_JOBID
echo $SLURM_NODELIST

# make sure conda on my path here
# (this is the version installed in my folder so we can get it on ssrde)
export PATH="/cube/neurocube/local/serenceslab/maggie/anaconda3/bin:$PATH"

echo $PATH
# activate my conda env here
# need to use full path, this is path to env on SSRDE
source activate /cube/neurocube/local/serenceslab/maggie/conda_envs/shape_dim

# change this path
ROOT=/cube/neurocube/local/serenceslab/maggie/shapeDim/

# put the code directory on your python path
PYTHONPATH=:${ROOT}Analysis/

# go to folder where script is located
cd ${ROOT}Analysis/

debug=0
n_threads=32
n_iter=1000;
# n_iter=10
rndseed=132334

python3 -c 'from binary_decoding import decode_binary; decode_binary.decode_withintask_permutationtest('${debug}', '${n_threads}', '${n_iter}', '${rndseed}',)'

# python3 -c 'from binary_decoding import decode_binary; decode_binary.decode_trainrepeat('${debug}', '${n_threads}')'
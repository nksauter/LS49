#!/bin/bash -l

#SBATCH -q special    # regular or special queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 00:15:00   # wall clock time limit
#SBATCH -J test_gpu_job
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH -G 1          # devices per node
#SBATCH -c 10         # total threads requested per node
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --exclusive

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

export USE_EXASCALE_API=True # "True" or "False" use granular host/device memory transfer
export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export N_SIM=30 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export ADD_BACKGROUND_ALGORITHM=cuda # cuda or jh or sort_stable
export CACHE_FHKL_ON_GPU=True # "True" or "False" use single object per rank
export DEVICES_PER_NODE=1
mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd;ls
srun -n 5 -c 2 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/step5_batch.py
echo "jobend $(date)";pwd;ls


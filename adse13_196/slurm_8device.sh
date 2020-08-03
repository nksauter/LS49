#!/bin/bash -l

#SBATCH -q special    # regular or special queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 01:10:00   # wall clock time limit
#SBATCH -J test_gpu_job
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH -G 8          # devices per node
#SBATCH -c 80         # total threads requested per node
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --exclusive

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

export LOG_BY_RANK=1 # Use Aaron's profiler/rank logger
export N_SIM=240 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export DEVICES_PER_NODE=8
mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
date;ls
srun -n 40 -c 2 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_161/step5_batch.py
date;ls
sleep 5

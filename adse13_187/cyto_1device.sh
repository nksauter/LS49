#!/bin/bash -l

#SBATCH -q special    # regular or special queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 02:05:00   # wall clock time limit
#SBATCH -J test_gpu_job
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH -G 1          # devices per node
#SBATCH -c 10         # total threads requested per node
#SBATCH -o job%j.out
#SBATCH -e job%j.err
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --exclusive

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export ADD_BACKGROUND_ALGORITHM=cuda # cuda or jh or sort_stable
export CACHE_FHKL_ON_GPU=True # "True" or "False" use single object per rank

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd;ls
srun -n 1 -c 2 libtbx.python ../cyto_batch.py N_total=1
echo "jobend $(date)";pwd;ls


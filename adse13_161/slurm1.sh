#!/bin/bash -l

#SBATCH -q regular    # premium queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 00:10:00   # wall clock time limit
#SBATCH -J test_gpu_job
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH --gres=gpu:8  # devices per node
#SBATCH -c 80         # total threads requested per node
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export LS49_BIG_DATA=/global/cscratch1/sd/nksauter/proj-h0918/ls49_big_data
export N_SIM=800 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export DEVICES_PER_NODE=8

srun -n 40 -c 2 libtbx.python ../modules/LS49/adse13_161/step5_batch.py

#!/bin/bash -l
#SBATCH -N 4                # Number of nodes
#SBATCH -J test_frame
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A m3890_g          # allocation
#SBATCH -C gpu
#SBATCH -q early_science    # regular queue
#SBATCH -t 01:00:00         # wall clock time limit
#SBATCH -n 64               # number of tasks
##SBATCH --ntasks-per-node=4
#SBATCH -c 8               # cpus per task
#SBATCH --ntasks-per-gpu=4
##SBATCH --gpus-per-task=1
##SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -o job%j.out
#SBATCH -e job%j.err

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
srun libtbx.python $MODULES/LS49/adse13_196/test_mpi.py

echo "jobend $(date)";pwd

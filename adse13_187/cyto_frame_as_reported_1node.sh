#!/bin/bash -l
#SBATCH -q special    # regular or special queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 00:15:00   # wall clock time limit
#SBATCH -J test_frame
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH -G 8          # devices total (not per node)
#SBATCH -c 80         # total threads requested per node
#SBATCH -o job%j.out
#SBATCH -e job%j.err
#SBATCH --exclusive
# test progression: 1 device+1 thread; 1 device+5 or 10 thread; 8 device

export WORK=$SCRATCH/adse13_187/bleededge/work
export BERNINA=$SCRATCH/adse13_187/bernina # location of data files
export MODULES=$SCRATCH/adse13_187/bleededge/alcc-recipes/cctbx/modules
export CFSX=/global/cfs/cdirs/m3562/swissfel
export CCTBX_DEVICE_PER_NODE=1
export DIFFBRAGG_USE_CUDA=1
export N_START=0
cd $WORK

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
#srun -n 40 -c 5 libtbx.python $MODULES/LS49/adse13_196/test_mpi.py
# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

# MINIMUM n >= 2 to implement the client-server pattern
srun -n 16 -c 5 --ntasks-per-gpu=2 --cpu-bind=cores libtbx.python $MODULES/LS49/adse13_187/cyto_frame_as_reported.py N_total=31
echo "jobend $(date)";pwd

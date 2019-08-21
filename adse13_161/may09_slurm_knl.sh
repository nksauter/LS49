#!/bin/bash -l

#SBATCH -q premium    # premium queue
#SBATCH -N 484        # Number of nodes
#SBATCH -t 12:00:00   # wall clock time limit
#SBATCH -J knl_job
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C knl,quad,cache
#SBATCH -A m2859      # allocation
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL

# starts off in /global/cscratch1/sd/nksauter/proj-h0918/test/052

# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 8 tasks per node for Haswell)

export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export LS49_BIG_DATA=/global/cscratch1/sd/nksauter/proj-h0918/ls49_big_data


srun -n 8228 -c 16 libtbx.python ../modules/LS49/sim/step5_batch.py

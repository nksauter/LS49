#!/bin/bash -l
#SBATCH -N 10               # Number of nodes
#SBATCH -J stills_proc
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A m3890_g          # allocation
#SBATCH -C gpu
#SBATCH -q early_science    # regular queue
#SBATCH -t 00:20:00         # wall clock time limit
#SBATCH -o job%j.out
#SBATCH -e job%j.err

export WORK=$SCRATCH/adse13_249/LY99
cd $WORK

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
#srun -n 32 -c 2 libtbx.python $MODULES/LS49/adse13_196/test_mpi.py
# one hour, 10 nodes should suffice
# NO PSF:
srun -n 320 -c 2 dials.stills_process ../index2.phil input.glob=../917302/LY99_MPIbatch_*.img.gz

# WITH PSF:
#srun -n 320 -c 2 dials.stills_process ../index2.phil input.glob=../917162/LY99_MPIbatch_*.img.gz
echo "jobend $(date)";pwd

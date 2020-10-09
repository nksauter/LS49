#!/bin/bash
#BSUB -P CHM137
#BSUB -W 00:20
#BSUB -nnodes 307
#BSUB -alloc_flags "gpumps nvme"
#BSUB -o job%J.out
#BSUB -e job%J.err
export WORK=/gpfs/alpine/chm137/scratch/nksauter/adse13_196
cd $WORK
mkdir ${LSB_JOBID}
cd ${LSB_JOBID}
export BBPATH=/mnt/bb/$USER/

jsrun -n 1842 -a 6 -c 6 -r 6 -g 1 python program.py
# Runs on 307 nodes, 1842 resource groups, 11052 MPI ranks
# Produces 90 GB in 122104 files collectively on the 307 burst buffers
# program only takes 138 seconds

jsrun -n 307 -r 1 find ${BBPATH} -type f -exec mv {} . \;
# Attempt to copy all the BB files back to GPFS takes 282 seconds
# Would like to cut this time down to a few seconds.  How?


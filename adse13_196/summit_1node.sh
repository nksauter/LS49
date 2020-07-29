#!/bin/bash
#BSUB -P CHM137
#BSUB -W 02:00
#BSUB -nnodes 1
#BSUB -alloc_flags gpumps

#BSUB -o job%j.out
#BSUB -e job%j.err

export N_SIM=240 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export DEVICES_PER_NODE=6
mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
date;ls
jsrun -n 1 -a 18 -g 6 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_161/step5_batch.py
date;ls

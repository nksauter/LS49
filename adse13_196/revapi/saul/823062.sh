#!/bin/bash -l
#SBATCH -N 10            # Number of nodes
#SBATCH -J NESAP_bench
#SBATCH -L SCRATCH       # job requires SCRATCH files
#SBATCH -A m3890_g       # allocation
#SBATCH -C gpu
#SBATCH -q early_science # regular or special queue
#SBATCH -t 00:07:00      # wall clock time limit
#SBATCH --gpus-per-node=4
#SBATCH -o job%j.out
#SBATCH -e job%j.err

export WORK=$SCRATCH/adse13_249/nesap
export BERNINA=$SCRATCH/adse13_187/bernina # location of data files
export MODULES=$SCRATCH/adse13_249/alcc-recipes/cctbx/modules
export CFSX=${SCRATCH}/cfs/cdirs/m3562/swissfel
export HOTPIXEL_MASK=$BERNINA/hopper_help_files/newmask_withbad.pkl
export STRUCTURE_FACTORS_MTZ_NAME=$BERNINA/hopper_help_files/100shuff.mtz
export CCTBX_DEVICE_PER_NODE=1
export DIFFBRAGG_USE_CUDA=1
export N_START=0
cd $WORK

export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export N_SIM=12800 # total number of images to simulate
export ADD_BACKGROUND_ALGORITHM=cuda # cuda or jh or sort_stable
export DEVICES_PER_NODE=4
export MOS_DOM=25

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
srun -n 160 -c 2 libtbx.python $MODULES/LS49/adse13_196/revapi/step5_batch.py
echo "jobend $(date)";pwd


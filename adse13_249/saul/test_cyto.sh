#!/bin/bash -l
#SBATCH -N 32               # Number of nodes
#SBATCH -J test_frame
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A m3890_g          # allocation
#SBATCH -C gpu
#SBATCH -q early_science    # regular queue
#SBATCH -t 01:00:00         # wall clock time limit
#SBATCH -n 128               # number of tasks
#SBATCH --ntasks-per-node=4
#SBATCH -c 32               # cpus per task
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -o job%j.out
#SBATCH -e job%j.err

export WORK=$SCRATCH/adse13_249/work
export BERNINA=$SCRATCH/adse13_187/bernina # location of data files
export MODULES=$SCRATCH/adse13_249/alcc-recipes/cctbx/modules
export CFSX=${SCRATCH}/cfs/cdirs/m3562/swissfel
export HOTPIXEL_MASK=$BERNINA/hopper_help_files/newmask_withbad.pkl
export STRUCTURE_FACTORS_MTZ_NAME=$BERNINA/hopper_help_files/100shuff.mtz
export CCTBX_DEVICE_PER_NODE=1
export DIFFBRAGG_USE_CUDA=1
export N_START=0
cd $WORK

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
srun libtbx.python $MODULES/LS49/adse13_196/test_mpi.py
srun libtbx.python $MODULES/LS49/adse13_187/cyto_frame_as_reported.py N_total=7534
echo "jobend $(date)";pwd

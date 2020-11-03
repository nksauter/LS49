#!/bin/bash
#BSUB -P CHM137
#BSUB -W 00:20
#BSUB -nnodes 115
#BSUB -alloc_flags "gpumps nvme"
#BSUB -o job%J.out
#BSUB -e job%J.err
cd $WORK
mkdir ${LSB_JOBID}
cd ${LSB_JOBID}

export BBPATH=/mnt/bb/$USER/
export USE_EXASCALE_API=True # "True" or "False" use granular host/device memory transfer
export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export N_SIM=99360 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export ADD_BACKGROUND_ALGORITHM=cuda # cuda or jh or sort_stable
export CACHE_FHKL_ON_GPU=True # "True" or "False" use single object per rank
export DEVICES_PER_NODE=1

echo "jobstart $(date)";pwd;ls
jsrun -n 690 -a 6 -c 6 -r 6 -g 1 libtbx.python \
  $(libtbx.find_in_repositories LS49)/adse13_196/step5_batch.py logger.outdir=${BBPATH}

#echo "verify burst buffer files"
#jsrun -n 115 -r 1 ls ${BBPATH}
#echo "verify local files"
#ls -l

jsrun -n 115 -r 1 find ${BBPATH} -type f -exec cp {} . \;

#echo "check that the files have been written to local GPFS"
#ls -l

echo "jobend $(date)";pwd
# -alloc_flags "gpumps" The GPU devices can be accessed by multiple MPI ranks
# DEVICES_PER_NODE=1    In the Summit architecture, the GPU DeviceID is always 0
# -n 6                  --nrs number of resource sets, 6 per node
# -a 7                  --tasks_per_rs number of MPI ranks per resource set
# -c 7                  --cpu_per_rs CPU cores per resource set (42 per node)
# -r 6                  --rs_per_host resource sets per node
# -d packed             --launch_distribution, packed is the default
# Explained at https://docs.olcf.ornl.gov/systems/summit_user_guide.html#batch-scripts

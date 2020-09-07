#!/bin/bash
#mkdir results
#cd results

export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export OMP_NUM_THREADS=60
export N_SIM=40 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=NKS # cuda or JH or NKS
export DEVICES_PER_NODE=1
echo "jobstart $(date)";pwd;ls
libtbx.python $(libtbx.find_in_repositories LS49)/adse13_161/step5_batch.py
echo "jobend $(date)";pwd


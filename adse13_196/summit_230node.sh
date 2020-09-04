#!/bin/bash
#BSUB -P CHM137
#BSUB -W 00:20
#BSUB -nnodes 230
#BSUB -alloc_flags "gpumps"
#BSUB -o job%J.out
#BSUB -e job%J.err
cd $WORK
mkdir ${LSB_JOBID}
cd ${LSB_JOBID}

export USE_EXASCALE_API=True # "True" or "False" use granular host/device memory transfer
export LOG_BY_RANK=1 # Use Aaron's profiler/rank logger
export N_SIM=100000 # total number of images to simulate
export ADD_SPOTS_ALGORITHM=cuda # cuda or JH or NKS
export DEVICES_PER_NODE=1
date;pwd;ls
#jsrun -n 6 -a 7 -c 7 -r 6 -g 1 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/test_mpi.py
jsrun -n 1380 -a 7 -c 7 -r 6 -g 1 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/step5_batch.py
#jsrun -n 6 -a 7 -c 7 -r 6 -g 1 libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/test_mpi.py
date;pwd;ls
# -alloc_flags "gpumps" The GPU devices can be accessed by multiple MPI ranks
# DEVICES_PER_NODE=1    In the Summit architecture, the GPU DeviceID is always 0
# -n 6                  --nrs number of resource sets, 6 per node
# -a 7                  --tasks_per_rs number of MPI ranks per resource set
# -c 7                  --cpu_per_rs CPU cores per resource set (42 per node)
# -r 6                  --rs_per_host resource sets per node
# -d packed             --launch_distribution, packed is the default
# Explained at https://docs.olcf.ornl.gov/systems/summit_user_guide.html#batch-scripts

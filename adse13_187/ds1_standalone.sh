#!/bin/bash -l
#SBATCH -q special    # regular or special queue
#SBATCH -N 1          # Number of nodes
#SBATCH -t 03:20:00   # wall clock time limit
#SBATCH -J test_frame
#SBATCH -L SCRATCH    # job requires SCRATCH files
#SBATCH -C gpu
#SBATCH -A m1759      # allocation
#SBATCH -G 8          # devices total (not per node)
#SBATCH -c 80         # total threads requested per node
#SBATCH -o job%j.out
#SBATCH -e job%j.err
# not #SBATCH --exclusive
# test progression: 1 device+1 thread; 1 device+5 or 10 thread; 8 device

export WORK=$SCRATCH/adse13_187/bleededge/work
export BERNINA=$SCRATCH/adse13_187/bernina # location of data files
export MODULES=$SCRATCH/adse13_187/bleededge/alcc-recipes/cctbx/modules
export CFSX=/global/cfs/cdirs/m3562/swissfel
export CCTBX_DEVICE_PER_NODE=1
export DIFFBRAGG_USE_CUDA=1
cd $WORK

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID
echo "jobstart $(date)";pwd
#srun -n 16 -c 5 libtbx.python $MODULES/LS49/adse13_196/test_mpi.py
# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 40 tasks per node for Cori-gpu)

# MINIMUM n >= 2 to implement the client-server pattern
#srun -n 8 -c 2 libtbx.python $MODULES/LS49/adse13_187/cyto_frame.py N_total=7
#srun -n 16 -c 5 --cpu-bind=cores --ntasks-per-gpu=2 libtbx.python $MODULES/LS49/adse13_187/cyto_frame.py N_total=7534
#srun -n 1 -c 5 libtbx.python $MODULES/LS49/adse13_187/cyto_frame.py N_total=1
srun -n 8 -c 5 --ntasks-per-gpu=1 --cpu-bind=cores libtbx.python $MODULES/LS49/adse13_187/cyto_frame.py N_total=750
#srun -n 8 -c 10 --cpu-bind=cores --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 libtbx.python $MODULES/LS49/adse13_187/cyto_frame.py N_total=15
echo "jobend $(date)";pwd

#tomorrow,
# figure out how to run the jobs without graphical output DONE
# scale up to 4 devices to figure out device ids. 6115 OK CCTBX delegated, LUNUS always device=0 (JBLOCKS=32)
# scale up to 8 devices, 10 min., 7 images. 6129. out of memory! I cannot stack more than one task per device
# scale up to 8 devices, 10 min., 7 images. use gpu binding. 6176 OK, configured every device as 0 (did not pass device id to lunus)
# scale up to 8 devices, 8 ranks for 1 hour, 70 images. 6188, 6202. Reveal some programming bugs
# scale up to 4 nodes for 4 hours to model 1280 images
#CUDA error 2 [/usr/common/software/sles15_cgpu/cuda/11.1.1/bin/../targets/x86_64-linux/include/cub/util_allocator.cuh, 446]: out of memory
#CUDA error 2 [/usr/common/software/sles15_cgpu/cuda/11.1.1/bin/../targets/x86_64-linux/include/cub/util_allocator.cuh, 490]: out of memory
#CUDA error 2 [/global/cscratch1/sd/nksauter/adse13_187/bleededge/alcc-recipes/cctbx/modules/lunus/lunus/cuda/lsort.cu, 109]: out of memory
# bugs dials.image_viewer down
# C++ bug, double-free
# Paley assertions fail
# lunus requires device-id==0
# lunus takes enormous memory; want to share the device
# monday morning can't even get any output at all.
# Fri:  make a separate execution for DS1.  Get before/after rmsd.  Implement tasks per device

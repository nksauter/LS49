#!/bin/bash -l
#SBATCH -N 32         #Use N nodes
#SBATCH -t 23:50:00  #Set wall clock time limit
#SBATCH -p regular  #Submit to the regular or debug or premium 'partition', or debug
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -A lcls
#SBATCH --job-name=step4K
#SBATCH --mail-user=nksauter@lbl.gov
#SBATCH --mail-type=ALL

# starts off in /global/cscratch1/sd/nksauter/proj-e/LS49

module load python/2.7-anaconda
source ../miniconda/bin/activate myEnv
source ../buildomp/setpaths.sh

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_DISPLAY_ENV=True
# -n, tasks to run; -N number of nodes; -c cpus per task;
# n = N x tasks_per_node (should be 8 tasks per node for Haswell)
srun -n 256 -c 8 libtbx.python ../modules/LS49/sim/step4Kbatch.py

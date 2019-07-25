#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 200
#SBATCH -t 04:00:00
#SBATCH -J my_job
#SBATCH -L SCRATCH
#SBATCH -C knl
#SBATCH -A lcls
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err

#srun -n 68 -c 2 libtbx.python ../modules/LS49/work2_for_aca_lsq/abc_background.py
export OMP_NUM_THREADS=8
export JSON_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export IMAGE_GLOB=/global/cscratch1/sd/nksauter/proj-h0918/HASWELL1/step6_MPIbatch_0%05d.img.gz
export USE_POSTREFINE=False
export MODEL_MODE=coarse_ground_truth #dials_refine #" | "coarse_ground_truth"

#srun -n 6800 -c 8
libtbx.python ../modules/LS49/work2_for_aca_lsq/abc_background.py

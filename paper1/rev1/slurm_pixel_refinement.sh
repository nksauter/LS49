#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 400
#SBATCH -t 48:00:00  # will be 04:00:00
#SBATCH -J my_job
#SBATCH -L SCRATCH
#SBATCH -C knl
#SBATCH -A lcls
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err

export WORK=/global/cscratch1/sd/nksauter/proj-paper1/work
export LS49_BIG_DATA=/global/cscratch1/sd/nksauter/proj-h0918/ls49_big_data
export OMP_NUM_THREADS=16
export BOOST_ADAPTBX_FPE_DEFAULT=1
#export NRANK=64
export JSON_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=dials_refine
export ABC_GLOB_A=/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_dials_refine/abcX%06d.pickle
export ABC_GLOB_C=/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_coarse_ground_truth/abcX%06d.pickle
export ABC_GLOB=/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_dials_refine/abcX%06d.pickle
export ABC_GLOB_PIXEL_REF=/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_pixel_refine/abcX%06d.pickle
export IMAGE_GLOB=/global/cscratch1/sd/nksauter/proj-h0918/HASWELL1/step6_MPIbatch_%06d.img.gz

# actually edit nanoBragg_nks to set num threads to 16.

# number of nodes (200) * 272 (knl) / c-stride(16) = number of ranks (3400)
srun -n 13600 -c 8 libtbx.python ../modules/LS49/ML_push/pixel_refinement.py


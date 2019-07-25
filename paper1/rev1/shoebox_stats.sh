#!/bin/bash -f

ln -s ../abc_coverage_pixel_refine ./abc_coverage

export JSON_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=dials_refine

srun -n 34 -c 8 libtbx.python ../modules/LS49/ML_push/shoebox_stats.py

rm ./abc_coverage


export WORK=/dev/shm/sauter
export LS49_BIG_DATA=${WORK}/ls49_big_data
export OMP_NUM_THREADS=64
export BOOST_ADAPTBX_FPE_DEFAULT=1
export NRANK=64
export JSON_GLOB=/net/dials/raid1/sauter/LS49_step6/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/net/dials/raid1/sauter/LS49_step6/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=dials_refine
export ABC_GLOB=/net/dials/raid1/sauter/paper1/abc_coverage_dials_refine/abcX%06d.pickle
export IMAGE_GLOB=/net/dials/raid1/sauter/LS49_step6/HASWELL1/step6_MPIbatch_%06d.img.gz
libtbx.python ../modules/LS49/ML_push/pixel_refinement.py

#!/bin/bash -l

export JSON_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=pixel_refine
export ABC_GLOB=/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_pixel_refine/abcX%06d.pickle

# negative control
libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=False \
LLG_evaluator.max_calls=11 \
LLG_evaluator.spoilHKL=True \
LLG_evaluator.title=ox_red_spoilHKL_50ka cohort=0

# red-red
libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_reduced_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=red_red_50ka

#!/bin/bash -l
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J my_job
#SBATCH -L SCRATCH
#SBATCH -C knl
#SBATCH -A m2859
#SBATCH -o slurm_replot_method4.out
#SBATCH -e slurm_replot_method4.err

export JSON_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=pixel_refine
export ABC_GLOB=/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_pixel_refine/abcX%06d.pickle

# ground truth, job 1
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot_Figure9.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_red_50k


# metal metal
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot_Figure9.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_50ka

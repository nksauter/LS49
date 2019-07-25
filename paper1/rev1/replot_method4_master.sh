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
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_red_50k

# negative control
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=False \
LLG_evaluator.max_calls=11 \
LLG_evaluator.spoilHKL=True \
LLG_evaluator.title=ox_red_spoilHKL_50ka cohort=0

# red-red
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_reduced_model starting_model.preset.FE2=Fe_reduced_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=red_red_50ka

# red-ox, job 2
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_reduced_model starting_model.preset.FE2=Fe_oxidized_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=red_ox_50k

# ox-ox
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_oxidized_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_ox_50ka

# metal metal
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_50ka

# fewer images, job 3
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=25000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_25k

srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=12000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_12k

srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=6000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_6k

# job 4
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=3000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_3k

srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=1500 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_1500

srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=750 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_750

# disjoint
srun -n 1 -c 16 libtbx.python ../modules/LS49/ML_push/replot.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=50000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
cohort=1 \
LLG_evaluator.title=metal_metal_c1_50k

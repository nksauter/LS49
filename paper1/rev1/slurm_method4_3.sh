#!/bin/bash -l
#SBATCH -q regular
#SBATCH -N 16
#SBATCH -t 48:00:00
#SBATCH -J my_job
#SBATCH -L SCRATCH
#SBATCH -C knl
#SBATCH -A m2859
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
export NRANK=64
export JSON_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/method3/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=pixel_refine
export ABC_GLOB=/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_pixel_refine/abcX%06d.pickle
# number of nodes (16) * 272 (knl) / c-stride(8) = number of ranks (1360)

# fewer images, job 3
srun -n 544 -c 8 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=25000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_25k > metal_metal_25k.out 2>&1

srun -n 544 -c 8 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=12000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_12k > metal_metal_12k.out 2>&1

srun -n 544 -c 8 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=6000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_6k > metal_metal_6k.out 2>&1

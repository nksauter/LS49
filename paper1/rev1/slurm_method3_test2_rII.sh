#!/bin/bash -l
#SBATCH -q premium
#SBATCH -N 80
#SBATCH -t 08:00:00
#SBATCH -J my_job
#SBATCH -L SCRATCH
#SBATCH -C knl
#SBATCH -A lcls
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
export NRANK=64
export JSON_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/global/cscratch1/sd/nksauter/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=dials_refine
export ABC_GLOB=/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_dials_refine/abcX%06d.pickle
# number of nodes (40) * 272 (knl) / c-stride(8) = number of ranks (1360)
# ground truth
srun -n 2720 -c 8 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=32000 \
LLG_evaluator.restraints_II_enable=True \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_red_32k_rII > ox_red_32k_rII.out 2>&1

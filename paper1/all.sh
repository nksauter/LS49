#!/bin/bash -f

export NRANK=64
# restraints or no restraints:
export JSON_GLOB=/net/dials/raid1/sauter/LS49_step6/integration/idx-step6_MPIbatch_0%05d.img_integrated_experiments.json
export PICKLE_GLOB=/net/dials/raid1/sauter/LS49_step6/integration/idx-step6_MPIbatch_0%05d.img_integrated.pickle
export USE_POSTREFINE=False
export MODEL_MODE=dials_refine
export ABC_GLOB=/net/dials/raid1/sauter/paper1/abc_coverage_coarse_ground_truth/abcX%06d.pickle

# negative control
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=False \
LLG_evaluator.max_calls=11 \
LLG_evaluator.spoilHKL=True \
LLG_evaluator.title=ox_red_spoilHKL_32k cohort=0 > ox_red_spoilHKL_32k.out 2>&1

# ground truth
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_reduced_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_red_32k > ox_red_32k.out 2>&1

# red-red
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_reduced_model starting_model.preset.FE2=Fe_reduced_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=red_red_32k > red_red_32k.out 2>&1

# red-ox
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_reduced_model starting_model.preset.FE2=Fe_oxidized_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=red_ox_32k > red_ox_32k.out 2>&1

# ox-ox
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_oxidized_model starting_model.preset.FE2=Fe_oxidized_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=ox_ox_32k > ox_ox_32k.out 2>&1

# metal metal
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_32k > metal_metal_32k.out 2>&1

# fewer images
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=16000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_16k > metal_metal_16k.out 2>&1

mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=8000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_8k > metal_metal_8k.out 2>&1

mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=4000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_4k > metal_metal_4k.out 2>&1

mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=2000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_2k > metal_metal_2k.out 2>&1

mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=1000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_1k > metal_metal_1k.out 2>&1

mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=500 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
LLG_evaluator.title=metal_metal_500 > metal_metal_500.out 2>&1

# disjoint
mpirun -c ${NRANK} libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py \
starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=32000 \
LLG_evaluator.enable_plot=False LLG_evaluator.plot_interpolation=True \
LLG_evaluator.max_calls=24 \
cohort=1 \
LLG_evaluator.title=metal_metal_c1_32k > metal_metal_c1_32k.out 2>&1

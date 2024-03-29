#!/bin/bash -l
#SBATCH -N 4                # Number of nodes
#SBATCH -J roi
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A m3890_g          # allocation
#SBATCH -C gpu
#SBATCH -q early_science    # regular queue
#SBATCH -t 00:10:00         # wall clock time limit
#SBATCH -o job%j.out
#SBATCH -e job%j.err

export WORK=$SCRATCH/adse13_249/LY99
cd $WORK

mkdir -p $SLURM_JOB_ID; cd $SLURM_JOB_ID

export TRIAL=ly99sim
export OUT_DIR=${PWD}
# NO PSF:
export DIALS_OUTPUT=${WORK}/927185
# WITH PSF:
#export DIALS_OUTPUT=${WORK}/927187

export DIFFBRAGG_USE_BRAGG=1

echo "dispatch.step_list = input balance filter statistics_unitcell model_statistics annulus
input.path=${DIALS_OUTPUT}
input.experiments_suffix=.img_integrated.expt
input.reflections_suffix=.img_integrated.refl
input.keep_imagesets=True
input.persistent_refl_cols=shoebox
input.persistent_refl_cols=xyzcal.px
input.persistent_refl_cols=xyzobs.px.value
input.persistent_refl_cols=delpsical.rad
input.persistent_refl_cols=panel
input.parallel_file_load.method=uniform
scaling.model=${WORK}/1m2a.pdb
scaling.unit_cell=67.2 59.8 47.2 90 110.3 90
scaling.space_group=C2
scaling.resolution_scalar=0.993420862158964
filter.algorithm=unit_cell
filter.unit_cell.algorithm=cluster
filter.unit_cell.cluster.covariance.file=${WORK}/covariance_ly99sim_30000.pickle
filter.unit_cell.cluster.covariance.component=0
filter.unit_cell.cluster.covariance.mahalanobis=4.0
filter.outlier.min_corr=-1.0
merging.d_max=None
merging.d_min=2.1
statistics.annulus.d_max=2.5
statistics.annulus.d_min=2.1
spread_roi.enable=True
spread_roi.strong=2.0
output.output_dir=${OUT_DIR}/annulus
diffBragg.fix.detz_shift=True
diffBragg.fix.eta_abc=False
diffBragg.logging.other_ranks_level=high
diffBragg.spectrum_from_imageset = True
diffBragg.downsamp_spec.skip=True
diffBragg.simulator.structure_factors.mtz_name=${WORK}/928123/out/ly99sim_all.mtz
diffBragg.simulator.structure_factors.mtz_column=IMEAN,SIGIMEAN
diffBragg.space_group=C2
diffBragg.method=L-BFGS-B
diffBragg.simulator.oversample=1
diffBragg.refiner.adu_per_photon=1.0
diffBragg.use_restraints=True
" > annulus.phil
# how can I log the minimizer step by step in real time.  It is a complete black box.
mkdir -p ${OUT_DIR}/${TRIAL}/out
mkdir -p ${OUT_DIR}/${TRIAL}/tmp

echo "jobstart $(date)";pwd
srun -n 128 -c 2 cctbx.xfel.merge annulus.phil
echo "jobend $(date)";pwd

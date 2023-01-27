#!/bin/bash -l
#SBATCH -N 1               # Number of nodes
#SBATCH -J annulus1
#SBATCH -L SCRATCH         # job requires SCRATCH files
#SBATCH -A m3562          # allocation
#SBATCH -C cpu
#SBATCH -q regular          # regular queue
#SBATCH -t 00:10:00         # wall clock time limit
#SBATCH -o %j.out
#SBATCH -e %j.err
mkdir -p $SLURM_JOB_ID; cd $SLURM_JOB_ID
export CCTBX_NO_UUID=1
export DIFFBRAGG_USE_CUDA=1
export CUDA_LAUNCH_BLOCKING=1
export NUMEXPR_MAX_THREADS=128
export SLURM_CPU_BIND=cores # critical to force ranks onto different cores.
#                             verify with ps -o psr <pid>
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export SIT_PSDM_DATA=/global/cfs/cdirs/lcls/psdm-sauter
#export SIT_PSDM_DATA=/pscratch/sd/p/psdatmgr/data/pmscr
export CCTBX_GPUS_PER_NODE=1
export XFEL_CUSTOM_WORKER_PATH=$MODULES/psii_spread/merging/application

echo "
input {
  path=/global/cfs/cdirs/m3562/users/vidyagan/p20231/common/results/r0080/000_rg002/out
  path=/global/cfs/cdirs/m3562/users/vidyagan/p20231/common/results/r0081/000_rg002/out
  experiments_suffix=_refined.expt
  reflections_suffix=_indexed.refl
  parallel_file_load.method=uniform
  parallel_file_load.balance=global1
}
dispatch.step_list = input balance filter statistics_unitcell model_statistics annulus
input.keep_imagesets=True
input.read_image_headers=False
input.persistent_refl_cols=shoebox
input.persistent_refl_cols=bbox
input.persistent_refl_cols=xyzcal.px
input.persistent_refl_cols=xyzcal.mm
input.persistent_refl_cols=xyzobs.px.value
input.persistent_refl_cols=xyzobs.mm.value
input.persistent_refl_cols=xyzobs.mm.variance
input.persistent_refl_cols=delpsical.rad
input.persistent_refl_cols=panel
input.parallel_file_load.method=uniform

filter.algorithm=unit_cell
filter.unit_cell.algorithm=cluster
filter.unit_cell.cluster.covariance.file=$CFS/lcls/sauter/LY99/p20231/uc_filter/covariance_run080_cells.pickle
filter.unit_cell.cluster.covariance.component=0
filter.unit_cell.cluster.covariance.mahalanobis=2.0

scaling.model=/global/cfs/cdirs/m3562/references/5tis.pdb
scaling.unit_cell=117.72 223.23 310.50 90.00  90.00  90.00
scaling.space_group=P212121
scaling.resolution_scalar=0.993420862158964

merging.d_min=3.0

statistics.annulus.d_max=3.5
statistics.annulus.d_min=3.0

spread_roi.enable=True
spread_roi.strong=1.0
output.log_level=0
exafel.trusted_mask=/sf/bernina/data/p20231/work/common/mask/69plus.mask
exafel.scenario=1

output.output_dir=out
output.save_experiments_and_reflections=True
" > trial8.phil
echo "jobstart $(date)";pwd
srun -n 10 -c 4 cctbx.xfel.merge trial8.phil
echo "jobend $(date)";pwd


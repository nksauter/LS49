cd /global/cfs/cdirs/lcls/sauter/LY99//ready2/3724210/out

export OUTPUT_FOLDER=1

mkdir -p $OUTPUT_FOLDER; cd $OUTPUT_FOLDER

export CCTBX_NO_UUID=1
export DIFFBRAGG_USE_CUDA=1
export CUDA_LAUNCH_BLOCKING=1
export NUMEXPR_MAX_THREADS=128
export SLURM_CPU_BIND=cores # critical to force ranks onto different cores. verify with ps -o psr <pid>
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export SIT_DATA=/global/common/software/lcls/psdm/data
export SIT_ROOT=/reg/g/psdm
export SIT_PSDM_DATA=/global/cfs/cdirs/lcls/psdm-sauter
export CCTBX_GPUS_PER_NODE=1
export XFEL_CUSTOM_WORKER_PATH=$MODULES/psii_spread/merging/application # User must export $MODULES path
export WERK=/global/cfs/cdirs/lcls/sauter/LY99/

echo "
dispatch.step_list = input annulus
input.path=$WERK/ready2/3724210/out
input.experiments_suffix=trial8_scenario3A_000000.expt
input.reflections_suffix=trial8_scenario3A_000000.refl
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

scaling.model=/global/cfs/cdirs/m3562/users/vidyagan/6ydi.pdb
scaling.unit_cell=107.00  107.00  304.01  90.00  90.00  90.00
scaling.space_group=P41212
scaling.resolution_scalar=0.993420862158964
scaling.pdb.k_sol=0.435

filter.unit_cell.cluster.covariance.file=$WERK/reference/covariance_run145_cells.pickle
filter.unit_cell.cluster.covariance.component=0

merging.d_max=None
merging.d_min=2.5

statistics.annulus.d_max=4.0
statistics.annulus.d_min=2.5

spread_roi.enable=True
# spread_roi.strong=1.0 # only use for initial annulus definition, not subsequent

output.log_level=1 # 0 = stdout stderr, 1 = terminal
output.output_dir=out
output.prefix=trial8_scenario3A
output.save_experiments_and_reflections=True

exafel.scenario=S1
exafel.static_fcalcs.path=$WERK/reference/mmo_static_fcalcs.pickle
exafel.static_fcalcs.whole_path=$WERK/reference/mmo_miller_array.pickle
exafel.static_fcalcs.action=read
exafel.trusted_mask=$WERK/reference/epix.mask
exafel.shoebox_border=0
exafel.context=kokkos_gpu
exafel.model.plot=False
exafel.model.mosaic_spread.value=0.0512
exafel.model.Nabc.value=48,48,24
exafel.debug.lastfiles=False # write out *.h5, *.mask for each image
exafel.debug.verbose=False
exafel.debug.finite_diff=-1
exafel.debug.eps=1.e-8
exafel.debug.format_offset=0
exafel.debug.energy_offset_eV=0
exafel.debug.energy_stride_eV=4.00
exafel.skin=False # whether to use diffBragg
exafel{
  refpar{
    label = *background *G
    background {
      algorithm=rossmann_2d_linear
      scope=spot
      slice_init=border
      slice=all
    }
    G {
      scope=lattice
      reparameterize=bound
    }
  }
}
exafel.metal=MMO2
sauter20.LLG_evaluator.enable_plot=True
sauter20.LLG_evaluator.title=tell
sauter20.LLG_evaluator.restraints.fp.sigma=0.04
sauter20.LLG_evaluator.restraints.fdp.sigma=0.08
sauter20.LLG_evaluator.max_calls=30
trumpet.plot_all.enable=False
trumpet.plot_all.savepng=True
" > trial8.phil
echo "jobstart $(date)";pwd
cctbx.xfel.merge trial8.phil
echo "jobend $(date)";pwd

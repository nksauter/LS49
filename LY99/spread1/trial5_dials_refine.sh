#!/bin/bash -l
#SBATCH -N 24               # Number of nodes on Perlmutter
#SBATCH -J test_spread
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A lcls             # allocation
#SBATCH -C cpu
#SBATCH -q regular          # regular queue
#SBATCH -t 00:20:00         # wall clock time limit
#SBATCH -o %j.out
#SBATCH -e %j.err

mkdir -p $SLURM_JOB_ID; cd $SLURM_JOB_ID

export CCTBX_NO_UUID=1
export DIFFBRAGG_USE_CUDA=1
export CUDA_LAUNCH_BLOCKING=1
export NUMEXPR_MAX_THREADS=128
export SLURM_CPU_BIND=cores # critical to force ranks onto different cores. verify with ps -o psr <pid>
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export SIT_PSDM_DATA=/global/cfs/cdirs/lcls/psdm-sauter
export CCTBX_GPUS_PER_NODE=1
export XFEL_CUSTOM_WORKER_PATH=$MODULES/psii_spread/merging/application # User must export $MODULES path
export WERK=/global/cfs/cdirs/lcls/sauter/LY99/

echo "
dispatch.step_list = input balance annulus trumpet
input.path=${WERK}/all_plots_bugreport/out
input.experiments_suffix=.expt
input.reflections_suffix=.refl
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

scaling.model=$WERK/reference/6ydi.pdb
scaling.unit_cell=107.00  107.00  304.01  90.00  90.00  90.00
scaling.space_group=P41212
scaling.resolution_scalar=0.993420862158964

filter.unit_cell.cluster.covariance.file=$WERK/reference/covariance_run145_cells.pickle
filter.unit_cell.cluster.covariance.component=0

merging.d_max=None
merging.d_min=2.1

statistics.annulus.d_max=2.9
statistics.annulus.d_min=2.5

spread_roi.enable=True
spread_roi.strong=1.0

output.log_level=0 # 0 = stdout stderr, 1 = terminal
output.output_dir=out
output.prefix=trial5_production
output.save_experiments_and_reflections=True

exafel.scenario=1
trumpet.plot_all.enable=True
trumpet.plot_all.savepng=True
trumpet.dials_refine=True
trumpet.debug_rank_0=False
trumpet.spectrum.recalibration=1.95
trumpet.spectrum.nbins=60
trumpet.spectrum.range=7060.0,7180.0
trumpet.residuals.cutoff=5.
trumpet.outlier.stddev=3.
trumpet.outlier.deff=2000.
trumpet.outlier.ucell=2.
trumpet.outlier.tff=0
trumpet.outlier.remove=True
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector {
      fix_list = Tau1,Tau2,Tau3
      hierarchy_level = 0
    }
    crystal {
      # allow refinement
      unit_cell {
        restraints {
          tie_to_target {
            values = 106.93,106.93,303.72,90,90,90
            sigmas = 0.04, 0.04, 0.15, 0.0, 0.0, 0.0 # use 1/2 std dev from the covariance file
            id = None # apply to all
          }
        }
      }
    }
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
" > trial5.phil
echo "jobstart $(date)";pwd
srun -n 96 -c 4 cctbx.xfel.merge trial5.phil
echo "jobend $(date)";pwd

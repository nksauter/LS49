export MODULES=$SCRATCH/dials/20220303/alcc-recipes/cctbx/modules
export CCTBX_DEVICE_PER_NODE=4
source $MODULES/../activate.sh

export LOG_BY_RANK=1 # Use Aaron's rank logger
export RANK_PROFILE=0 # 0 or 1 Use cProfiler, default 1
export N_SIM=100 # total number of images to simulate
export ADD_BACKGROUND_ALGORITHM=cuda # always cuda for any GPU
export DEVICES_PER_NODE=4
export MOS_DOM=25
rm -r data; mkdir data; pushd data
echo "jobstart $(date)";pwd
srun -n 32 -G 4 libtbx.python $MODULES/LS49/adse13_196/revapi/LY99_batch.py context=kokkos_gpu noise=True psf=False attenuation=True oversample=1
echo "jobend $(date)";pwd
popd


echo "output {
  composite_output = False
}
dispatch {
  index=True
  refine=True
  integrate=True
}
mp {
  method = mpi
}
geometry {
  detector {
    panel {
      origin = -169.048 169.049 -141.7
    }
  }
}
spotfinder {
  filter {
    min_spot_size = 3
  }
  threshold {
    dispersion {
      gain = 1.0 # for nanoBragg sim
      sigma_background=2
      sigma_strong=2
      global_threshold=10
      kernel_size=6 6
    }
  }
  filter.d_min=1.9 # gives 10085=60 spots
}
indexing {
  stills.refine_candidates_with_known_symmetry=True
  known_symmetry {
    space_group = C2
    unit_cell = 67.2 59.8 47.2 90 110.2 90
  }
}
integration {
  background.simple.outlier.plane.n_sigma=10
  debug.output=True
  debug.separate_files=False
  summation {
    detector_gain = 1.0 # for nanoBragg sim
  }
}
profile.gaussian_rs.centroid_definition=com
#dispatch.pre_import=True # would otherwise import every nanoBragg image in every rank
#output.composite_output=True
output.logging_dir=. # demangle by rank
" > dials.phil

rm -r dials_out; mkdir dials_out; rm -r dials_log; mkdir dials_log
srun -n 32 dials.stills_process data/*.img.gz dials.phil output.output_dir=dials_out output.logging_dir=dials_log

source /pscratch/sd/d/dwpaley/dials/20220303/alcc-recipes/cctbx/activate.sh #FIXME
rm -r out

PROCESSED_DIR=$PWD/dials_out

export TRIAL=ly99sim
export CCTBX_NO_UUID=1

echo "dispatch.step_list = input filter statistics_unitcell model_statistics annulus
input.path=${PROCESSED_DIR}
input.experiments_suffix=.img_integrated.expt
input.reflections_suffix=.img_integrated.refl
input.keep_imagesets=True
input.read_image_headers=False
input.persistent_refl_cols=shoebox
input.persistent_refl_cols=bbox
input.persistent_refl_cols=xyzcal.px
input.persistent_refl_cols=xyzobs.px.value
input.persistent_refl_cols=delpsical.rad
input.persistent_refl_cols=panel
input.parallel_file_load.method=uniform
output.output_dir=out
scaling.model=/pscratch/sd/d/dwpaley/ly99/pdb/1m2a.pdb
scaling.unit_cell=67.2 59.8 47.2 90 110.3 90
scaling.space_group=C2
scaling.resolution_scalar=0.993420862158964
filter.algorithm=unit_cell
filter.unit_cell.algorithm=cluster
filter.unit_cell.cluster.covariance.file=/pscratch/sd/d/dwpaley/ly99/cov/covariance_1304534.cells.pickle
filter.unit_cell.cluster.covariance.component=0
filter.unit_cell.cluster.covariance.mahalanobis=4.0
filter.outlier.min_corr=-1.0
merging.d_max=None
merging.d_min=2.1
statistics.annulus.d_max=2.5
statistics.annulus.d_min=2.1
spread_roi.enable=True
spread_roi.strong=2.0
output.log_level=0 # stdout stderr
exafel.trusted_mask=/pscratch/sd/d/dwpaley/ly99/mask/pixels.mask
exafel.scenario=ds1
exafel.shoebox_border=0
exafel.context=kokkos_gpu
exafel.model.plot=False
exafel.model.Nabc.value=50,50,50
exafel.debug.lastfiles=False
exafel.debug.verbose=False


diffBragg {

  logging.rank0_level=high
  logging.other_ranks_level=high

  no_Nabc_scale=False
  method="L-BFGS-B"
  use_restraints=False
  space_group=C2
  spectrum_from_imageset = True
  downsamp_spec {
    skip = True
  }
  roi {
    shoebox_size=12
    fit_tilt=True
    fit_tilt_using_weights = False
    hotpixel_mask = None
    reject_edge_reflections = False
    reject_roi_with_hotpix = False
    pad_shoebox_for_background_estimation=10
    mask_outside_trusted_range=True
  }
  refiner {
    adu_per_photon = 1
    sigma_r=3
  }
  simulator {
    oversample=4
    crystal.has_isotropic_ncells = False
    structure_factors.from_pdb.name=/pscratch/sd/d/dwpaley/ly99/pdb/1m2a.pdb
    structure_factors.from_pdb.k_sol=0.435
    structure_factors.from_pdb.b_sol=46.0
    init_scale = 1
    beam.size_mm = 0.001
    detector.force_zero_thickness = True
  }
  init {
    Nabc=[72,72,72]
    G=100
  }
  mins {
    Nabc=[3,3,3]
    detz_shift=-1.5
    RotXYZ=[-15,-15,-15]
    G=0
  }
  maxs {
    RotXYZ=[15,15,15]
    Nabc=[1600,1600,1600]
    G=1e12
    detz_shift=1.5
  }
  sigmas {
    RotXYZ=[1e-3,1e-3,1e-3]
  }
  fix.detz_shift=True
  refiner.num_devices=4
  logging.log_refined_params=True
}

" > annulus.phil

echo "jobstart $(date)";pwd
DIFFBRAGG_USE_CUDA=1 srun -n 32 -G 4 cctbx.xfel.merge annulus.phil
echo "jobend $(date)";pwd

from __future__ import division, print_function

# DiffBragg stage 1, report version 2 dated Oct 4, 2021
ds1_params_v2="""no_Nabc_scale=False
temp=0.01
stepsize=0.05
ucell_edge_perc=15
ucell_ang_abs=1
niter=0
niter_per_J=1
rescale_params=True
quiet=True
space_group=P6522
use_restraints = False
method="L-BFGS-B"
spectrum_from_imageset = True
downsamp_spec {
  delta_en = 1
}
roi {
  shoebox_size = 13
  fit_tilt=True
  fit_tilt_using_weights = False
  hotpixel_mask = /global/cscratch1/sd/nksauter/adse13_187/bleededge/work/hopper_help_files/newmask_withbad.pkl
  reject_edge_reflections = False
  reject_roi_with_hotpix = False
  pad_shoebox_for_background_estimation=10
}
refiner {
  adu_per_photon = 9.481
  sigma_r=1.5
}
simulator {
  total_flux=1e12
  oversample=2
  crystal.has_isotropic_ncells = False
  structure_factors.mtz_name = /global/cscratch1/sd/nksauter/adse13_187/bleededge/work/hopper_help_files/100shuff.mtz
  structure_factors.mtz_column = "F(+),F(-)"
  init_scale = 1
  beam.size_mm = 0.001
  detector.force_zero_thickness = True
}
init {
  Nabc=[31.03,30.41,10.63]
  G=100
  B=10
}
mins {
  Nabc=[3,3,3]
  detz_shift=-1.5
  B=0
  RotXYZ=[-15,-15,-15]
  G=0
}
maxs {
  RotXYZ=[15,15,15]
  Nabc=[1600,1600,1600]
  G=1e8
  detz_shift=1.5
  B=1e5
}
sigmas {
  RotXYZ=[1e-3,1e-3,1e-3]
}
fix.detz_shift=True
outdir="."
logging.other_ranks_level="high"
"""

# DiffBragg stage 1, report version 4 dated Oct 6, 2021
# improved radial offset, improved correlation between delta R and delta Psi
import os
ds1_params_v4="""
ucell_edge_perc=15
ucell_ang_abs=1
method="Nelder-Mead"
#method="L-BFGS-B"
use_diffuse_models=True
#nelder_mead_maxfev=None
#nelder_mead_fatol=1000
#niter=20
spectrum_from_imageset = True
downsamp_spec {
  delta_en = 0.25
}
roi {
  fit_tilt=True
  fit_tilt_using_weights = False
  # supply mask file on command line
  hotpixel_mask = %(HOTPIXEL_MASK)s
  reject_edge_reflections = False
  pad_shoebox_for_background_estimation=10
}
refiner {
  adu_per_photon = 9.481
  sigma_r=10
}
simulator {
  oversample=4
  structure_factors.mtz_name = %(STRUCTURE_FACTORS_MTZ_NAME)s
  structure_factors.mtz_column = "F(+),F(-)"
  beam.size_mm = 0.001
  detector.force_zero_thickness = True
}
init {
  Nabc = 50.000000 50.000000 37.500000
  G = 10.000000
  diffuse_sigma = 1 1 1
  #diffuse_sigma = 0.583037 0.458631 0.704636
  diffuse_gamma = 100 100 100
  #diffuse_gamma = 238.873078 168.346065 73.587935
}
mins {
  detz_shift=-1.5
  RotXYZ=[-15,-15,-15]
}
maxs {
  detz_shift = 1.5
  Nabc = 1600 1600 1600
  RotXYZ = 15 15 15
  G = 100000
}
sigmas {
  RotXYZ=[1e-3,1e-3,1e-3]
}
use_restraints = True
betas {
  detz_shift = 1e-08
  ucell = 0.001 0.001
  RotXYZ = 1e-05
  Nabc = 50.000000 50.000000 50.000000
  G = 1000.000000
}
centers {
  ucell = 78.61 265.12
  Nabc = 50.000000 50.000000 37.500000
  G = 10.000000
}
fix.detz_shift=True
fix.diffuse_gamma=False
fix.diffuse_sigma=False
outdir="."
"""%(os.environ)


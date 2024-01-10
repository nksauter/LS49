from __future__ import division

import os

from libtbx import test_utils
import libtbx.load_env

tst_list = [
  "$D/tests/tst_numpy_lsq.py",
  "$D/tests/tst_spectrum_iterator.py",
  "$D/tests/tst_structure_factors.py",
  "$D/tests/tst_sf_energies.py",
  "$D/tests/tst_mosaic_orientations.py",
  "$D/tests/tst_crystal_orientations.py",
  "$D/tests/tst_jh_add_spots.py",
]

OpenMP_optional = [
    "$D/tests/tst_monochromatic_image.py",  #OpenMP (optional)
    ["$D/tests/tst_cpu_add_spots_double_precision.py","mono"],  #OpenMP (optional)
  ]

OpenMP_required = [
    "$D/tests/tst_polychromatic_image.py",  #OpenMP (required)
    ["$D/tests/tst_cpu_add_spots_double_precision.py","poly"],  #OpenMP (required)
  ]
tst_list_parallel = []
tst_list_parallel_slow = []
if libtbx.env.build_options.enable_openmp_if_possible:
  tst_list_parallel += OpenMP_optional
  tst_list_parallel_slow += OpenMP_required
else:
  tst_list += OpenMP_optional

prepend = []
OPT = libtbx.env.build_options
   # these three tests, break portability after realizing that the spectral dispersion curve
   # comes directly from the nexus master and the expt must be read with check_format=True

# check if /global exists before adding tests that depend on files in /global
if os.path.isdir('/global'):
  if OPT.enable_cuda:  prepend.append(
    ["$D/adse13_187/cyto_batch.py", "N_total=1", "test_pixel_congruency=True",
      "mosaic_method=double_random",
      "mosaic_spread_samples=50", "write_output=False", "test_without_mpi=True",
      "log.outdir=mp1c",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=cuda"
    ])
  if OPT.enable_kokkos:  prepend.append(
    ["$D/adse13_187/cyto_batch.py", "N_total=1", "test_pixel_congruency=True",
      "mosaic_method=double_random",
      "mosaic_spread_samples=50", "write_output=False", "test_without_mpi=True",
      "log.outdir=mp1k",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=kokkos_gpu"
    ])
  if OPT.enable_cuda:  prepend.append(
    ["$D/adse13_187/tst_multipanel_argchk.py", "N_total=1",
      "mosaic_spread_samples=50", "test_without_mpi=True", "log.outdir=mp2c",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=cuda"
    ])
  if OPT.enable_kokkos:  prepend.append(
    ["$D/adse13_187/tst_multipanel_argchk.py", "N_total=1",
      "mosaic_spread_samples=50", "test_without_mpi=True", "log.outdir=mp2k",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=kokkos_gpu"
    ])
  if OPT.enable_cuda:  prepend.append(
    ["$D/adse13_187/tst_write_file_action.py", "N_total=1",
      "mosaic_method=double_random",
      "mosaic_spread_samples=50", "test_without_mpi=True", "log.outdir=mp3c",
      "write_output=False",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=cuda"
    ])
  if OPT.enable_kokkos:  prepend.append(
    ["$D/adse13_187/tst_write_file_action.py", "N_total=1",
      "mosaic_method=double_random",
      "mosaic_spread_samples=50", "test_without_mpi=True", "log.outdir=mp3k",
      "write_output=False",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "context=kokkos_gpu"
    ])
  # this test pair is developmental only, not portable, only checks run-without-crash, not results:
  if OPT.enable_cuda:  prepend.append(
    ["$D/adse13_187/tst_write_file_action.py", "N_total=1", "write_output=True", "write_experimental_data=True",
      "mosaic_spread_samples=62", "test_without_mpi=True", "log.outdir=mp4c",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "mask_file=/global/cfs/cdirs/m3562/nks/adse13_187/13_221/event_648.mask",
      "context=cuda"
    ])
  if OPT.enable_kokkos:  prepend.append(
    ["$D/adse13_187/tst_write_file_action.py", "N_total=1", "write_output=True", "write_experimental_data=True",
      "mosaic_spread_samples=62", "test_without_mpi=True", "log.outdir=mp4k",
      "nxmx_local_data=/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5",
      "mask_file=/global/cfs/cdirs/m3562/nks/adse13_187/13_221/event_648.mask",
      "context=kokkos_gpu"
    ])
# end tests that depend on files in /global

if OPT.enable_cuda:
  prepend = prepend + [
   "$D/adse13_196/tst_gpu_channels.py",
   "$D/adse13_196/revapi/tst_step5_batch_single_process_GPU.py",
  ]
tst_list_parallel = prepend + tst_list_parallel
if  OPT.enable_cuda:
  tst_list_parallel = tst_list_parallel + [
  ["$D/tests/tst_cuda_add_spots.py","mono"],
  ["$D/tests/tst_cuda_add_spots.py","poly"],
  # Laue hasn't worked recently, return to the issue later ["$D/tests/tst_cuda_add_spots.py","laue"]
  ]

def run_standalones():
  build_dir = libtbx.env.under_build("LS49")
  dist_dir = libtbx.env.dist_path("LS49")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run_standalones()

""" Discussion.
tst_numpy_lsq.py:
  Shows that numpy least squares fit works properly.
  Fits the relationship between spectral peak and xtc-recorded energy.
  Numpy code is not thread safe; this workaround temporily sets OMP_NUM_THREADS to 1
tst_spectrum_iterator.py:
  sameness of processed spectra to reference run
tst_input_pdb.py:
  fetches PDB code 1m2a and checks identity to old reference file
tst_structure_factors.py:
  the computed structure factors without energy dependence, in both P1 and C2
tst_sf_energies.py:
  the computed structure factors at selected energies
tst_mosaic_orientations.py:
  the mosaic domains
tst_crystal_orientations.py
  the 100000 random orientations
tst_monochromatic_image.py
  monochromatic "quick" simulation of raw image 0
  the air and water scatterers
tst_polychromatic_image.py
  the final 100-channel polychromatic image
  includes the energy-dependent effects of Fe site absorption, for oxidized & reduced Fe
  includes separate verification of each-wavelength energy contribution to raw image
tst_jh_add_spots.py
  Test the original JH code before refactoring for OpenMP
  Equivalence of add_nanobragg_spots(original) and add_nanobragg_spots_nks(refactored)
Future:
break down Fmodel into bulk solvent, site-scatterer, anomalous site-scatterer, and redox-dependent anomalous site-scatterer
Write tests that replicate ACA-2018 Bayesian estimates.
"""

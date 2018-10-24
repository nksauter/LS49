from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/tests/tst_numpy_lsq.py",
  "$D/tests/tst_spectrum_iterator.py",
  "$D/tests/tst_input_pdb.py",
  "$D/tests/tst_structure_factors.py",
  "$D/tests/tst_sf_energies.py",
  "$D/tests/tst_mosaic_orientations.py",
  "$D/tests/tst_crystal_orientations.py",
  "$D/tests/tst_monochromatic_image.py",
  "$D/tests/tst_polychromatic_image.py",
  "$D/tests/tst_jh_add_spots.py",
  )

OPT = libtbx.env.build_options
if OPT.enable_cuda:
  tst_list = tst_list + (
  ["$D/tests/tst_cuda_add_spots.py","mono"],
  ["$D/tests/tst_cuda_add_spots.py","poly"],
  ["$D/tests/tst_cuda_add_spots.py","laue"])

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

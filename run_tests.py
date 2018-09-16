from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/tests/tst_numpy_lsq.py",
  "$D/tests/tst_spectrum_iterator.py",
  "$D/tests/tst_input_pdb.py",
  "$D/tests/tst_structure_factors.py",
  )

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
the computed structure factors at selected energies
the mosaic domains
the 100000 random orientations
the standard-wavelength raw image
a sample wavelength energy contribution to raw image
the air and water scatters
the final 100-channel image
Make sure that I can test the original JH code as well!!!
In other words, equivalence of add_nanobragg_spots() and add_nanobragg_spots_nks
()
Or alternatively, re-express my step5 problem in terms of his code.
Write tests that replicate my Bayesian estimates.
"""

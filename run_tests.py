from __future__ import division
from libtbx import test_utils
import libtbx.load_env

tst_list = (
  "$D/tests/tst_numpy_lsq.py",
  "$D/tests/tst_spectrum_iterator.py",
  "$D/tests/tst_input_pdb.py",
  )

def run_standalones():
  build_dir = libtbx.env.under_build("LS49")
  dist_dir = libtbx.env.dist_path("LS49")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run_standalones()

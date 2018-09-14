from __future__ import division, print_function
from six.moves import cPickle
import os
from LS49.spectra import generate_spectra
from libtbx.test_utils import approx_equal

def get_results():
  ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  R = cPickle.load(open(os.path.join(ls49_big_data,"data/spectra209.pickle"),"rb"))
  return R
generate_spectra.get_results = get_results

def run():
  from LS49.spectra.generate_spectra import linear_fit # calls numpy linalg lsqsq
  R = get_results()
  LF = linear_fit(R)
  reason = """Numpy fails to give correct slope and intercept if OMP_NUM_THREADS > 1.
  Test asserts that subject LS49 code provides a workaround to temporarily set OMP_NUM_THREADS to 1.
  Test assert that linear fit agrees with pre-computed values.
  Fit describes relationship between header energy and spectrometer peak."""
  assert approx_equal(LF.m, 0.0590617409556), reason
  assert approx_equal(LF.c, 7028.05604768), reason

if __name__=="__main__":
  run()
  print("OK")

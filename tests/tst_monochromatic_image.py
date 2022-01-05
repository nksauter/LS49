from __future__ import division, print_function
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
import six

from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data

def run_monochromatic():
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=True)

# temporarily increasing tolerances. Reset to delta=50, count=10
def compare_two_images(reference, test, tolerance_delta=2540, tolerance_count=3984027, verbose_pixels=False, verbose=True):
  print ("Comparing",reference,test)
  try:
    from dxtbx.format.Registry import Registry
    get_format = Registry.find
  except ImportError:
    from dxtbx.format import Registry
    get_format = Registry.get_format_class_for_file

  beam=[]
  data = []
  detector =[]
  headers = []
  for i,fk in enumerate([reference,test]):
    format_instance = get_format(fk)
    instance = format_instance(fk)
    beam.append( instance.get_beam() )
    detector.append( instance.get_detector() )
    data.append( instance.get_raw_data() )
    if verbose: #optional test
      print (beam[-1])
      print (instance.get_goniometer())
      print (detector[-1])
      print (instance.get_scan())
    headers.append( instance.get_smv_header(fk) )

  if  headers[0]==headers[1]:
    print ("Both headers identical")
  else:
    #print headers[0]
    #print headers[1]
    for key in headers[0][1]:
      if key not in headers[1][1]:  print ("second data lacks key",key)
      elif headers[0][1][key] != headers[1][1][key]:
        print ("Key comparison:",key,headers[0][1][key], headers[1][1][key])

  assert len(data[1])==len(data[0])
  diff_data = data[1]-data[0]
  from scitbx.array_family import flex
  abs_diff_data = flex.abs(diff_data)

  #verbose output of diffs
  if verbose_pixels:
    nonzero_deltas = (abs_diff_data > 0)
    ndiff = 0
    for idiff,diff in enumerate(diff_data):
      if diff!=0:
        if ndiff < 200: print ("difference index %d:(%d,%d)"%(idiff,idiff//3000,idiff%3000),diff,data[0][idiff]) # only print the first 200 differences
        ndiff += 1
    print("There are %d differences"%ndiff)

  # first filter, do not allow any |delta| greater than tolerance_delta
  N_large_deltas = (abs_diff_data > tolerance_delta).count(True)
  assert N_large_deltas == 0, "%d pixels have |delta| larger than cutoff %d with max %d"%(N_large_deltas,tolerance_delta, max(abs_diff_data))

  # next filter, allow a maximum |delta| of tolerance_delta for no more than tolerance_count pixels
  nonzero_deltas = (abs_diff_data > 0)
  N_non_zero_deltas = nonzero_deltas.count(True)
  if N_non_zero_deltas > 0:
    print("%d pixels have |delta| up to %d"%(N_non_zero_deltas,tolerance_delta))
  assert N_non_zero_deltas <= tolerance_count, "%d pixels have |delta| up to %d"%(N_non_zero_deltas,tolerance_delta)

# temporarily increasing tol. reset to 1E-7
def compare_two_raw_images(reference, test, tol=0.003): # TODO: run more tests to decide on the default tolerance
  from six.moves import cPickle as pickle
  from scitbx.array_family import flex
  with open(reference,'rb') as F:
    if six.PY3:
      ref_array = pickle.load(F, encoding="bytes")
    else:
      ref_array = pickle.load(F)
  with open(test,'rb') as F:
    test_array = pickle.load(F)
  print("\nComparing raw image: '%s' with the reference: '%s'"%(test, reference))
  diff_array = test_array - ref_array
  if diff_array.all_eq(0.0):
    print ("There are 0 differences\n")
  else:
    stats = flex.mean_and_variance(diff_array.as_1d()) # flex.mean_and_variance works only on 1d arrays
    diff_mean = stats.mean()
    diff_std = stats.unweighted_sample_standard_deviation()
    diff_min = flex.min(diff_array)
    diff_max = flex.max(diff_array)
    print("Differences: range (%.2E to %.2E); mean %.2E; std %.2E"%(diff_min, diff_max, diff_mean, diff_std))
    # assert acceptable differences
    assert abs(diff_mean) < tol, "The raw image is different from the reference with abs(diff_mean) %f > tol %f." % (diff_mean, tol)

if __name__=="__main__":
  run_monochromatic()
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000.img.gz"), test="./step5_000000.img.gz")
  # test the raw photons due to Bragg scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_001.img"), test="./step5_000000_intimage_001.img")
  # add in the effects of water scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_002.img"), test="./step5_000000_intimage_002.img")
  # add in the effects of air scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_003.img"), test="./step5_000000_intimage_003.img")
  print("OK")

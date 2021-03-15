from __future__ import division, print_function
from libtbx.development.timers import Profiler

#provision for reading masks
def mask_from_file(mask_file):
  P = Profiler("fast mask from file as boolean array")
  from scitbx.array_family import flex
  import pickle
  with open(mask_file,"rb") as M:
    mask = pickle.load(M)
  monolithic_mask = flex.bool()
  for ipnl in range(len(mask)):
    pnl = mask[ipnl]
    monolithic_mask.extend(pnl.as_1d())
  assert len(monolithic_mask)==256*254*254
  return monolithic_mask

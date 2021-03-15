from __future__ import division, print_function

#provision for reading masks
def mask_from_file(mask_file):
  from libtbx.development.timers import Profiler
  P = Profiler("mask from file as boolean array")
  from scitbx.array_family import flex
  import pickle
  with open(mask_file,"rb") as M:
    mask = pickle.load(M)
    monolithic_mask = flex.bool()
  for ipnl in range(len(mask)):
    pnl = mask[ipnl]
    for item in pnl: monolithic_mask.append(item)
  assert len(monolithic_mask)==256*254*254
  return monolithic_mask

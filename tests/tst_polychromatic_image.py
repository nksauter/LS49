from __future__ import division, print_function
import os
from six.moves import cPickle
from LS49.sim import debug_utils

class channel_extractor:
  def __init__(self):
    self.data=dict()
  def extract(self, channel_no, data):
    if channel_no%10 == 0:
      self.data[channel_no]=data
      print ("Storing channel information for energy %d"%channel_no)
debug_utils.channel_extractor = channel_extractor
filename="energy_dependent_diffraction.pickle"

from LS49.tests import tst_monochromatic_image
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment

def run_polychromatic(create):
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=False) #actually perform the 100-channel simulation
  from LS49.sim.step5_pad import CHDBG_singleton
  if create:
    cPickle.dump(CHDBG_singleton.data,
      open(os.path.join(ls49_big_data,"reference",filename),"wb"),cPickle.HIGHEST_PROTOCOL)
  else:
    CH_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference",filename),"rb"))
    assert len(CHDBG_singleton.data)==10,"Should have recorded 10 energy channels"
    for key in CHDBG_singleton.data:
        assert CHDBG_singleton.data[key] == CH_ref[key],"Energy-channel results should agree with reference"

if __name__=="__main__":
  run_polychromatic(create=False)
  tst_monochromatic_image.compare_two_images(
    reference=os.path.join(ls49_big_data,"reference","step5_MPIbatch_000000.img.gz"), test="./step5poly_000000.img.gz")
  print("OK")

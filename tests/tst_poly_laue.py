from __future__ import division, print_function
import os

from LS49.tests import tst_monochromatic_image
from LS49 import ls49_big_data

def run_polychromatic(create):
  #from LS49.sim import step5_laue
  #step5_laue.add_spots_algorithm = "JH"
  from LS49.sim.step5_laue import tst_all
  tst_all(quick=False) #actually perform the 100-channel simulation

if __name__=="__main__":
  run_polychromatic(create=False)
  tst_monochromatic_image.compare_two_images(
    reference=os.path.join(ls49_big_data,"reference","step5laue_remote_ref_000000.img.gz"),
    test="./step5laue_000000.img.gz",
    tolerance_count=200)
  print("OK")

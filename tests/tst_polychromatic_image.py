from __future__ import division, print_function
import os
from LS49.tests import tst_monochromatic_image
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment

def run_polychromatic():
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=False)

if __name__=="__main__":
  run_polychromatic()
  tst_monochromatic_image.compare_two_images(
    reference=os.path.join(ls49_big_data,"reference","step5_MPIbatch_000000.img.gz"), test="./step5poly_000000.img.gz")
  print("OK")

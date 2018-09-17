from __future__ import division, print_function
import os
from LS49.tests.tst_monochromatic_image import compare_two_images
from LS49.sim import step5_pad

step5_pad.add_spots_algorithm = "JH"
# key idea:  this swaps in the production (James Holton) add_nanoBragg_spots() implementation
# rather than the Nick Sauter specialization refactored for OpenMP.

def run_monochromatic():
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=True,prefix="jh_step5")

if __name__=="__main__":
  run_monochromatic()
  ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000.img.gz"), test="./jh_step5_000000.img.gz")
  # test the raw photons due to Bragg scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_001.img"), test="./jh_step5_000000_intimage_001.img")
  # add in the effects of water scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_002.img"), test="./jh_step5_000000_intimage_002.img")
  # add in the effects of air scatter:
  compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_003.img"), test="./jh_step5_000000_intimage_003.img")
  print("OK")

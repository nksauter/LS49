from __future__ import division, print_function
import os
from LS49.tests.tst_monochromatic_image import compare_two_images, compare_two_raw_images
from LS49.sim import step5_pad

step5_pad.add_spots_algorithm = "cuda"
# key idea:  this swaps in the production (James Holton) add_nanoBragg_spots() implementation
# rather than the Nick Sauter specialization refactored for OpenMP.

def run_monochromatic():
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=True,prefix="cuda_step5")

def run_polychromatic(create):
  from LS49.sim.step5_pad import tst_all
  tst_all(quick=False,prefix="cuda_step5") #actually perform the 100-channel simulation

def run_laue(create):
  from LS49.sim import step5_laue
  step5_laue.add_spots_algorithm = "cuda"
  from LS49.sim.step5_laue import tst_all
  tst_all(quick=False,prefix="cuda_step5") #perform the 100-channel simulation with single kernel call

#explanation:  for mono and poly, use the step5_pad simulation with separate structure factors for each energy channel
# for laue, we want a single kernel call, so use step5_laue that sets only a single value for structure factors,
# and an appropriate reference image calculated that way on CPU.

if __name__=="__main__":
  import sys
  mode = sys.argv[1]
  assert mode in ["mono","poly","laue"]
  if mode == "laue":
    run_laue(create=False)
    ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5laue_remote_ref_000000.img.gz"), test="./cuda_step5laue_000000.img.gz",tolerance_count=200)
  elif mode == "poly":
    run_polychromatic(create=False)
    ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_MPIbatch_000000.img.gz"), test="./cuda_step5poly_000000.img.gz")
    # compare raw photon values (Bragg scatter only):
    compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","step5poly_cuda_ref_dblprec_000000.pickle"), test="./cuda_step5poly_000000_dblprec_001.pickle")

  else:
    run_monochromatic()
    ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000.img.gz"), test="./cuda_step5_000000.img.gz")
    # test the raw photons due to Bragg scatter:
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_001.img"), test="./cuda_step5_000000_intimage_001.img")
    # add in the effects of water scatter:
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_002.img"), test="./cuda_step5_000000_intimage_002.img")
    # add in the effects of air scatter:
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5_000000_intimage_003.img"), test="./cuda_step5_000000_intimage_003.img")
  print("OK")

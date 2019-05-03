from __future__ import division, print_function
import os
from LS49.tests.tst_monochromatic_image import compare_two_images


if __name__=="__main__":
  from LS49.sim import step5_pad_chan_sim
  import sys

  algo = sys.argv[1]
  assert( algo in ["cuda","NKS","JH"])
  step5_pad_chan_sim.ADD_SPOTS_ALGORITHM = algo

  print ("\n\n<><><><>\nUSING ALGO: %s\n><><><><>" % algo)

  step5_pad_chan_sim.tst_all(prefix="cuda_step5_pad_chan_sim_")
  ls49_big_data = os.environ["LS49_BIG_DATA"]

  compare_two_images(
    reference=os.path.join(
            ls49_big_data,
            "reference",
            "cuda_step5_pad_chan_sim_poly_000000.img.gz"),
    test="./cuda_step5_pad_chan_sim_poly_000000.img.gz")

  print("OK")

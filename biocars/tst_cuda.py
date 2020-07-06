from __future__ import division, print_function
import os
from LS49.tests.tst_monochromatic_image import compare_two_images

def run_monochromatic():
  pass

def run_polychromatic(create):
  pass

def run_laue(create):
  from LS49.biocars import laue
  laue.add_spots_algorithm = "cuda"
  from LS49.biocars.laue import tst_all
  tst_all(quick=False,prefix="biocars")

if __name__=="__main__":
  import sys
  mode = sys.argv[1]
  assert mode in ["mono","poly","laue"]
  if mode == "laue":
    run_laue(create=False)
    #from LS49 import ls49_big_data
    #compare_two_images(reference=os.path.join(ls49_big_data,"reference","step5laue_remote_ref_000000.img.gz"), test="./cuda_step5laue_000000.img.gz",tolerance_count=200)
  print("OK")

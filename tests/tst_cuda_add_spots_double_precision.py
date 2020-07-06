from __future__ import division, print_function
import os
import sys
from LS49.tests.tst_monochromatic_image import compare_two_images, compare_two_raw_images
from LS49.sim import step5_pad

step5_pad.add_spots_algorithm = "cuda"

# This test uses production (James Holton) add_nanoBragg_spots() implementation for cuda.
# The reference is add_nanoBragg_spots() implementation refactored by Nick Sauter for OpenMP.
# The reference was created by tst_cpu_add_spots_double_precision.py
# The test compares only the photons due to Bragg scatter, no air or water scatter effects.

def run_monochromatic():
  from LS49.sim.step5_pad import tst_all
  prefix = "cuda_step5"
  tst_all(quick=True,prefix=prefix,save_bragg=True) # quick=True means: perform a single wavelength simulation

def run_polychromatic():
  from LS49.sim.step5_pad import tst_all
  prefix = "cuda_step5"
  tst_all(quick=False,prefix=prefix,save_bragg=True) # quick=False means: perform multiple wavelength simulation

if __name__=="__main__":
  import LS49.utils.safe_to_write as s2w
  from LS49 import ls49_big_data
  mode = sys.argv[1]
  assert mode in ["mono","poly"]
  if mode == "poly":
    # make sure we don't overwrite existing images
    s2w.cwd_safe_to_write(["*cuda_step5poly_000000_intimage*.img","*cuda_step5poly_000000.img.gz", "*cuda_step5poly_000000_dblprec*.pickle"])
    # create images
    run_polychromatic()
    # compare to reference integer image, smv-formatted
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_intimage_001.img"), test="cuda_step5poly_000000_intimage_001.img")
    # compare to reference double precision image, numpy array pickle
    compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_dblprec_001.pickle"), test="./cuda_step5poly_000000_dblprec_001.pickle")
  else:
    # make sure we don't overwrite existing images
    s2w.cwd_safe_to_write(["*cuda_step5_000000_intimage*.img","*cuda_step5_000000.img.gz", "*cuda_step5_000000_dblprec*.pickle"])
    # create images
    run_monochromatic()
    # compare to reference integer image, smv-formatted
    compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_intimage_001.img"), test="cuda_step5_000000_intimage_001.img")
    # compare to reference double precision image, numpy array pickle
    compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_dblprec_001.pickle"), test="./cuda_step5_000000_dblprec_001.pickle")
  print("OK")

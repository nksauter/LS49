from __future__ import division, print_function
import os
from LS49.tests.tst_monochromatic_image import compare_two_images, compare_two_raw_images

# This test uses production (James Holton) add_nanoBragg_spots() implementation refactored by Nick Sauter for OpenMP.
# The test compares only the raw photons due to Bragg scatter, no air or water scattering effects.

def run_monochromatic(create):
  from LS49.sim.step5_pad import tst_all
  prefix = "cpu_step5"
  if create:
    prefix = "ref_" + prefix
  tst_all(quick=True,prefix=prefix,save_bragg=True) # quick=True means: perform a single wavelength simulation

def run_polychromatic(create):
  from LS49.sim.step5_pad import tst_all
  prefix = "cpu_step5"
  if create:
    prefix = "ref_" + prefix
  tst_all(quick=False,prefix=prefix,save_bragg=True) # quick=False means: perform multiple wavelength channel simulation

if __name__=="__main__":
  import sys
  mode = sys.argv[1]
  assert mode in ["mono","poly"]
  create = False # change this flag to True to create a new reference image
  ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  if mode == "poly":
    run_polychromatic(create=create)
    if not create:
      # compare to reference integer image, smp-formatted
      compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_intimage_001.img"), test="cpu_step5poly_000000_intimage_001.img")
      # compare to reference double precision image calculated on CPU
      compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_dblprec_001.pickle"), test="./cpu_step5poly_000000_dblprec_001.pickle")
  else:
    run_monochromatic(create=create)
    if not create:
      # compare to reference integer image, smp-formatted
      compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_intimage_001.img"), test="cpu_step5_000000_intimage_001.img")
      # compare to reference double precision image calculated on CPU
      compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_dblprec_001.pickle"), test="./cpu_step5_000000_dblprec_001.pickle")
  print("OK")

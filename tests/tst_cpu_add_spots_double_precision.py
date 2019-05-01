from __future__ import division, print_function
import os
import sys
from LS49.tests.tst_monochromatic_image import compare_two_images, compare_two_raw_images

# This test uses production (James Holton) add_nanoBragg_spots() implementation refactored by Nick Sauter for OpenMP.
# The test compares only the photons due to Bragg scatter, no air or water scatter effects.

create_ref = False # change this flag to True to create new reference images. Place all references under $LS49_BIG_DATA/references.

def run_monochromatic(create_ref):
  from LS49.sim.step5_pad import tst_all
  prefix = "cpu_step5"
  if create_ref:
    prefix = "ref_" + prefix
  tst_all(quick=True,prefix=prefix,save_bragg=True) # quick=True means: perform a single wavelength simulation

def run_polychromatic(create_ref):
  from LS49.sim.step5_pad import tst_all
  prefix = "cpu_step5"
  if create_ref:
    prefix = "ref_" + prefix
  tst_all(quick=False,prefix=prefix,save_bragg=True) # quick=False means: perform multiple wavelength simulation

if __name__=="__main__":
  # make sure this test doesn't overwrite some existing images
  import LS49.utils.safe_to_write as s2w
  s2w.cwd_safe_to_write(["*cpu_step5*_000000_intimage*.img","*cpu_step5*_000000.img.gz", "*cpu_step5*_000000_dblprec*.pickle"])

  mode = sys.argv[1]
  assert mode in ["mono","poly"]
  ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  if mode == "poly":
    run_polychromatic(create_ref=create_ref)
    if not create_ref:
      # compare to reference integer image, smv-formatted
      compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_intimage_001.img"), test="cpu_step5poly_000000_intimage_001.img")
      # compare to reference double precision image, numpy array pickle
      compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5poly_000000_dblprec_001.pickle"), test="./cpu_step5poly_000000_dblprec_001.pickle")
  else:
    run_monochromatic(create_ref=create_ref)
    if not create_ref:
      # compare to reference integer image, smv-formatted
      compare_two_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_intimage_001.img"), test="cpu_step5_000000_intimage_001.img")
      # compare to reference double precision image, numpy array pickle
      compare_two_raw_images(reference=os.path.join(ls49_big_data,"reference","ref_cpu_step5_000000_dblprec_001.pickle"), test="./cpu_step5_000000_dblprec_001.pickle")
  print("OK")

from __future__ import division
from LS49.tests.tst_monochromatic_image import compare_two_images
import os, sys

image_simulation_directory = sys.argv[1]
assert os.path.isdir(image_simulation_directory)
simulation_template = os.path.join(image_simulation_directory, "step5_MPIbatch_%06d.img.gz")

from LS49 import ls49_big_data
reference_directory = os.path.join(ls49_big_data,"adse13_196","reference")
reference_template = os.path.join(reference_directory, "step5_MPIbatch_%06d.img.gz")

for i in range(240):

  newfile = simulation_template%i
  oldfile = reference_template%i

  print(oldfile,newfile)
  import os
  if not os.path.isfile(oldfile):continue
  if not os.path.isfile(newfile):continue
  print("COMPARISON",i)
  compare_two_images(reference=oldfile, test=newfile, verbose=False)

print("OK")

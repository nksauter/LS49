from __future__ import division, print_function
import pickle
import os,shutil
os.environ["JSON_GLOB"]="null"
os.environ["PICKLE_GLOB"]="null"
os.environ["USE_POSTREFINE"]="null"
os.environ["MODEL_MODE"]="null"
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
from scitbx.array_family import flex
abc_glob_pixel_refine = "/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_pixel_refine/abcX%06d.pickle"
abc_glob_pixel_refine_90percent = "/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_pixel_refine_90percent/abcX%06d.pickle"
deltafast = flex.double()
deltaslow = flex.double()
acceptance_dict = {}
with open("LLGabstract.txt","r") as F:
  for line in F:
    tokens=line.split()
    key = int(tokens[2])
    value = float(tokens[19])
    if value<0.04:  acceptance_dict[key]="OK"

print("There are %d accepted keys"%len(acceptance_dict))
idone=0
for key in range(100000):
  print(key)
  if key not in acceptance_dict: continue
  try:
    pixel = abc_glob_pixel_refine%key
    pixel90 = abc_glob_pixel_refine_90percent%key
    assert os.path.isfile(pixel)
  except Exception:
    continue

  shutil.copyfile(src=pixel, dst=pixel90)
  idone+=1
  print ("copied %d"%(idone),pixel90)

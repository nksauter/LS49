from __future__ import division, print_function
import pickle
import os
os.environ["JSON_GLOB"]="null"
os.environ["PICKLE_GLOB"]="null"
os.environ["USE_POSTREFINE"]="null"
os.environ["MODEL_MODE"]="null"
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
from scitbx.array_family import flex
abc_glob_dials_refine = "/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_dials_refine/abcX%06d.pickle"
abc_glob_pixel_refine = "/global/cscratch1/sd/nksauter/proj-paper1/abc_coverage_pixel_refine/abcX%06d.pickle"
from scitbx.matrix import col
deltafast = flex.double()
deltaslow = flex.double()
acceptance_dict = {}
with open("LLGabstract.txt","r") as F:
  for line in F:
    tokens=line.split()
    key = int(tokens[2])
    value = float(tokens[19])
    if value<0.04:  acceptance_dict[key]="OK"

for key in range(2000):
  print(key)
  if key not in acceptance_dict: continue
  try:
    pixel = pickle.load(open(abc_glob_pixel_refine%key,"rb"))
    dials = pickle.load(open(abc_glob_dials_refine%key,"rb"))
  except IOError:
    continue

  assert len(pixel) == len(dials)


  for ispot in range(len(pixel)):
    dials_roi = dials[ispot].roi
    pixel_roi = pixel[ispot].roi
    focus = dials_roi.focus()

    S = flex.double(range(focus[0]))
    F = flex.double(range(focus[1]))
    # matrix giving the slow coordinate:
    cslow = S.matrix_outer_product(flex.double([1]*focus[1]))
    # matrix giving the fast coordinate:
    cfast = flex.double([1]*focus[0]).matrix_outer_product(F)

    sum_dials_roi = flex.sum(dials_roi)
    dials_expectation_value = col( ( flex.sum(cfast*dials_roi)/sum_dials_roi ,
                                     flex.sum(cslow*dials_roi)/sum_dials_roi ) )

    sum_pixel_roi = flex.sum(pixel_roi)
    pixel_expectation_value = col( ( flex.sum(cfast*pixel_roi)/sum_pixel_roi ,
                                     flex.sum(cslow*pixel_roi)/sum_pixel_roi ) )

    delta_position = pixel_expectation_value - dials_expectation_value
    deltafast.append(delta_position[0])
    deltaslow.append(delta_position[1])
    print (delta_position.elems)
    #from IPython import embed; embed()

from matplotlib import pyplot as plt
plt.plot(deltaslow,deltafast,"b,")
plt.axes().set_aspect("equal")
plt.show()

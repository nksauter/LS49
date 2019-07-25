from __future__ import division, print_function
import pickle
import os
os.environ["JSON_GLOB"]="null"
os.environ["PICKLE_GLOB"]="null"
os.environ["USE_POSTREFINE"]="null"
os.environ["MODEL_MODE"]="null"
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
from scitbx.array_family import flex
abc_glob_dials_refine = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_dials_refine/abcX%06d.pickle"
abc_glob_pixel_refine = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_pixel_refine/abcX%06d.pickle"
abc_glob_coarse_ground_truth = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_coarse_ground_truth/abcX%06d.pickle"
from scitbx.matrix import col
deltafast = flex.double()
deltaslow = flex.double()

for key in range(10000):
  print(key)
  try:
    pixel = pickle.load(open(abc_glob_pixel_refine%key,"rb"))
    dials = pickle.load(open(abc_glob_coarse_ground_truth%key,"rb"))
  except IOError:
    continue

  if not len(pixel) == len(dials):
    # only one image in 10000 fails
    continue


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

from matplotlib import pyplot as plt
#from IPython import embed; embed()
statss = flex.mean_and_variance(deltaslow)
print("Slow axis stats, mean, stddev=",statss.mean(), statss.unweighted_sample_standard_deviation(),"on N=",len(deltaslow))
statsf = flex.mean_and_variance(deltafast)
print("Fast axis stats, mean, stddev=",statsf.mean(), statsf.unweighted_sample_standard_deviation())

fig, ax = plt.subplots()

# the histogram of the data
n, bins, patches = ax.hist(deltaslow, bins=40, normed=0, histtype="step", range=(-2.0,2.0), color='blue')
n, bins, patches = ax.hist(deltafast, bins=40, normed=0, histtype="step", range=(-2.0,2.0), color='red')

ax.set_xlabel('spot shift (pixels)')
ax.set_ylabel('Count')
ax.set_title(r'Histogram of slow, fast pixel shift')
plt.show()
# commented code for scatter plot
plt.plot(deltaslow,deltafast,"b,")
plt.axes().set_aspect("equal")
plt.show()

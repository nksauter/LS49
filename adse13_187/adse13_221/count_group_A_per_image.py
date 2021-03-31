from __future__ import division
import os
from dials.array_family import flex
dirpath = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b"
from dxtbx.model.experiment_list import ExperimentListFactory

def from_files():
  for x in range(7534):
    if x%500==0: print(x)
    expt_file = "split_%04d.expt"%x
    expts = ExperimentListFactory.from_json_file(os.path.join(dirpath,expt_file),
                                              check_format=False)
    refl_file = "split_%04d.refl"%x
    refl_table = flex.reflection_table.from_file(os.path.join(dirpath,refl_file))

    yield refl_table,expts[0]

def run_job():
  count_spots_per_image = flex.double()
  for ref,exp in from_files():
    idx = ref["miller_index"]
    count_spots_per_image.append(len(idx))

  if True:
      from matplotlib import pyplot as plt
      plt.hist(count_spots_per_image, range=(0,700), bins=70, rwidth=0.8, color="orange", log=True)
      plt.title("Histogram of Bragg spots on each image")
      plt.xlabel("Number of Bragg spots")
      plt.ylabel("Number of images")
      plt.show()


if __name__=="__main__":
  run_job()


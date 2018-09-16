from __future__ import division, print_function
from six.moves import cPickle, range
from scitbx.array_family import flex
import os

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
filename = "crystal_orientations.pickle"

def model(create):
  N_crystal_orientations = 100000
  random_orientations = []
  mersenne_twister = flex.mersenne_twister(seed=0)
  for iteration in range(N_crystal_orientations):
    random_orientations.append( mersenne_twister.random_double_r3_rotation_matrix() )

  if create: # write the reference for the first time
      cPickle.dump(random_orientations,
      open(os.path.join(ls49_big_data,"reference",filename),"wb"),cPickle.HIGHEST_PROTOCOL)
  else: # read the reference and assert sameness to production run
      ori_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference",filename),"rb"))
      for x in range(len(random_orientations)):
        assert random_orientations[x] == ori_ref[x]

if __name__=="__main__":
  model(create=False)
  print("OK")

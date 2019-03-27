from __future__ import division, print_function
from six.moves import cPickle, range
from scitbx.array_family import flex
import scitbx
import os
import math
from scitbx.matrix import col

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
filename = "mosaic_domains.pickle"

def channel_wavelength_fmodel(create):
  N_mosaic_domains = 25
  mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev

  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  scitbx.random.set_random_seed(1234)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=mosaic_spread_deg * math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(N_mosaic_domains)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )

  if create: # write the reference for the first time
      cPickle.dump(UMAT_nm,
      open(os.path.join(ls49_big_data,"reference",filename),"wb"),cPickle.HIGHEST_PROTOCOL)
  else: # read the reference and assert sameness to production run
      UMAT_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference",filename),"rb"))
      for x in range(len(UMAT_nm)):
        print(x," ".join(
          ["%18.15f"%UMAT_ref[x][z] for z in range(9)]
        ))
        assert UMAT_nm[x] == UMAT_ref[x]

if __name__=="__main__":
  channel_wavelength_fmodel(create=False)
  print("OK")

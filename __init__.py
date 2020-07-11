from __future__ import division
import libtbx.load_env
srcdir = libtbx.env.dist_path("LS49")
import os.path
modules_dir = os.path.dirname(os.path.abspath(srcdir))
ls49_big_data = os.path.join(modules_dir,"ls49_big_data")

def legacy_random_orientations(N_total=100000):
  from scitbx.array_family import flex
  mt = flex.mersenne_twister_legacy_boost_1_63(seed=0)
  random_orientations = []
  for iteration in range(N_total):
    random_orientations.append( mt.random_double_r3_rotation_matrix() )
  return random_orientations


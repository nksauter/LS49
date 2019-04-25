from __future__ import division, print_function
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import sqr
import libtbx.load_env # possibly implicit
from cctbx import crystal
from time import time
from omptbx import omp_get_num_procs

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%

# Develop procedure for MPI control

def tst_one(image,spectra,crystal,random_orientation):
  iterator = spectra.generate_recast_renormalized_image(image=image,energy=7120.,total_flux=1e12)

  quick = False
  if quick: prefix_root="step5_batch_%06d"
  else: prefix_root="step5_MPIbatch_%06d"

  file_prefix = prefix_root%image
  rand_ori = sqr(random_orientation)
  from LS49.sim.step5_pad import run_sim2smv
  run_sim2smv(prefix = file_prefix,crystal = crystal,spectra=iterator,rotation=rand_ori,quick=quick,rank=rank)

if __name__=="__main__":
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  import os,omptbx
  workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
  omptbx.omp_set_num_threads(workaround_nt)
  N_total = 100000 # number of items to simulate
  N_stride = size # total number of worker tasks
  print("hello from rank %d of %d"%(rank,size),"with omp_threads=",omp_get_num_procs())
  import datetime
  start_elapse = time()
  if rank == 0:
    print("Rank 0 time", datetime.datetime.now())
    from LS49.spectra.generate_spectra import spectra_simulation
    from LS49.sim.step5_pad import microcrystal
    print("hello2 from rank %d of %d"%(rank,size))
    SS = spectra_simulation()
    C = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
    mt = flex.mersenne_twister(seed=0)
    random_orientations = []
    for iteration in range(N_total):
      random_orientations.append( mt.random_double_r3_rotation_matrix() )
    transmitted_info = dict(spectra = SS,
                            crystal = C,
                            random_orientations = random_orientations)
  else:
    transmitted_info = None
  transmitted_info = comm.bcast(transmitted_info, root = 0)
  comm.barrier()
  parcels = list(range(rank,N_total,N_stride))
  while len(parcels)>0:
    import random
    idx = random.choice(parcels)
    cache_time = time()
    print("idx------start-------->",idx,"rank",rank,time())
    tst_one(image=idx,spectra=transmitted_info["spectra"],
            crystal=transmitted_info["crystal"],random_orientation=transmitted_info["random_orientations"][idx])
    parcels.remove(idx)
    print("idx------finis-------->",idx,"rank",rank,time(),"elapsed",time()-cache_time)
  print("OK exiting rank",rank,"at",datetime.datetime.now(),"seconds elapsed",time()-start_elapse)

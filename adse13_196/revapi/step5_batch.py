from __future__ import division, print_function
from time import time
start_elapse = time()

from six.moves import range
from scitbx.matrix import sqr
import libtbx.load_env # possibly implicit
from cctbx import crystal
from omptbx import omp_get_num_procs
from scitbx.array_family import flex

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.adse13_196.revapi import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
from LS49.sim.util_fmodel import gen_fmodel
from LS49.adse13_196.revapi.step5_pad import data
from simtbx import get_exascale
# %%%%%%

# Develop procedure for MPI control

# first, sfall_channels DONE
# later, spectra to spectra iter
# data, why are we reading those files in all ranks?
# sfall_main not used?
# evaluate air + water as a singleton

def parse_input():
  from iotbx.phil import parse
  master_phil="""
    logger {
      outdir = .
        .type = path
        .help = Use "/mnt/bb/${USER}" for Summit NVME burst buffer
    }
    context = kokkos_gpu kokkos_cpu *cuda
      .type = choice
      .help = backend for parallel execution
      .help = Note, for now the assumption that default==cuda is baked in to the tests
      .help = specifically tst_step5_batch_single_process_GPU.py
  """
  phil_scope = parse(master_phil)
  # The script usage
  import libtbx.load_env # implicit import
  help_message = '''ADSE13-196.'''
  usage = ""
  '''Initialize the script.'''
  from dials.util.options import OptionParser
  # Create the parser
  parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        epilog=help_message)

  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options

def tst_one(image,spectra,crystal,random_orientation,sfall_channels,gpu_channels_singleton,rank,params):

  iterator = spectra.generate_recast_renormalized_image(image=image%100000,energy=7120.,total_flux=1e12)

  quick = False
  if quick: prefix_root="step5_batch_%06d"
  else: prefix_root="step5_MPIbatch_%06d"

  file_prefix = prefix_root%image
  rand_ori = sqr(random_orientation)
  from LS49.adse13_196.revapi.step5_pad import run_sim2smv
  run_sim2smv(prefix = file_prefix,
              crystal = crystal,
              spectra=iterator,rotation=rand_ori,quick=quick,rank=rank,
              gpu_channels_singleton=gpu_channels_singleton,
              sfall_channels=sfall_channels,params=params)

def run_step5_batch(test_without_mpi=False):
  params,options = parse_input()
  log_by_rank = bool(int(os.environ.get("LOG_BY_RANK",0)))
  rank_profile = bool(int(os.environ.get("RANK_PROFILE",1)))
  if log_by_rank:
    import io, sys
  if rank_profile:
    import cProfile
    pr = cProfile.Profile()
    pr.enable()

  if test_without_mpi:
    from LS49.adse13_196.mock_mpi import mpiEmulator
    MPI = mpiEmulator()
  else:
    from libtbx.mpi4py import MPI

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  import omptbx
  workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
  omptbx.omp_set_num_threads(workaround_nt)
  N_total = int(os.environ["N_SIM"]) # number of items to simulate
  N_stride = size # total number of worker tasks
  print("hello from rank %d of %d"%(rank,size),"with omp_threads=",omp_get_num_procs())
  import datetime
  start_comp = time()

  # now inside the Python imports, begin energy channel calculation

  wavelength_A = 1.74 # general ballpark X-ray wavelength in Angstroms
  wavlen = flex.double([12398.425/(7070.5 + w) for w in range(100)])
  direct_algo_res_limit = 1.7

  local_data = data() # later put this through broadcast

  GF = gen_fmodel(resolution=direct_algo_res_limit,
                  pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()

  # Generating sf for my wavelengths
  sfall_channels = {}
  for x in range(len(wavlen)):
    if rank > len(wavlen): break
    if x%size != rank: continue

    GF.reset_wavelength(wavlen[x])
    GF.reset_specific_at_wavelength(
                     label_has="FE1",tables=local_data.get("Fe_oxidized_model"),newvalue=wavlen[x])
    GF.reset_specific_at_wavelength(
                     label_has="FE2",tables=local_data.get("Fe_reduced_model"),newvalue=wavlen[x])
    sfall_channels[x]=GF.get_amplitudes()

  reports = comm.gather(sfall_channels, root = 0)
  if rank==0:
    sfall_channels = {}
    for report in reports:  sfall_channels.update(report)
  comm.barrier()

  print(rank, time(), "finished with the calculation of channels, now construct single broadcast")

  if rank == 0:
    print("Rank 0 time", datetime.datetime.now())
    from LS49.spectra.generate_spectra import spectra_simulation
    from LS49.adse13_196.revapi.step5_pad import microcrystal
    print("hello2 from rank %d of %d"%(rank,size))
    SS = spectra_simulation()
    C = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
    from LS49 import legacy_random_orientations
    random_orientations = legacy_random_orientations(N_total)
    transmitted_info = dict(spectra = SS,
                            crystal = C,
                            sfall_info = sfall_channels,
                            random_orientations = random_orientations)
  else:
    transmitted_info = None
  transmitted_info = comm.bcast(transmitted_info, root = 0)
  comm.barrier()
  parcels = list(range(rank,N_total,N_stride))

  print(rank, time(), "finished with single broadcast, now set up the rank logger")

  if log_by_rank:
    expand_dir = os.path.expandvars(params.logger.outdir)
    log_path = os.path.join(expand_dir,"rank_%d.log"%rank)
    error_path = os.path.join(expand_dir,"rank_%d.err"%rank)
    #print("Rank %d redirecting stdout/stderr to"%rank, log_path, error_path)
    sys.stdout = io.TextIOWrapper(open(log_path,'ab', 0), write_through=True)
    sys.stderr = io.TextIOWrapper(open(error_path,'ab', 0), write_through=True)

  print(rank, time(), "finished with the rank logger, now construct the GPU cache container")

  import random
  gpu_instance = get_exascale("gpu_instance", params.context)
  gpu_energy_channels = get_exascale("gpu_energy_channels", params.context)

  gpu_run = gpu_instance( deviceId = rank % int(os.environ.get("DEVICES_PER_NODE",1)) )

  gpu_channels_singleton = gpu_energy_channels (
    deviceId = gpu_run.get_deviceID())
    # singleton will instantiate, regardless of gpu, device count, or exascale API

  comm.barrier()
  while len(parcels)>0:
    idx = random.choice(parcels)
    cache_time = time()
    print("idx------start-------->",idx,"rank",rank,time())
    # if rank==0: os.system("nvidia-smi")
    tst_one(image=idx,spectra=transmitted_info["spectra"],
        crystal=transmitted_info["crystal"],
        random_orientation=transmitted_info["random_orientations"][idx],
        sfall_channels=transmitted_info["sfall_info"], gpu_channels_singleton=gpu_channels_singleton,
        rank=rank,params=params
    )
    parcels.remove(idx)
    print("idx------finis-------->",idx,"rank",rank,time(),"elapsed",time()-cache_time)
  comm.barrier()
  del gpu_channels_singleton
  # avoid Kokkos allocation "device_Fhkl" being deallocated after Kokkos::finalize was called
  print("Overall rank",rank,"at",datetime.datetime.now(),"seconds elapsed after srun startup %.3f"%(time()-start_elapse))
  print("Overall rank",rank,"at",datetime.datetime.now(),"seconds elapsed after Python imports %.3f"%(time()-start_comp))
  if rank_profile:
    pr.disable()
    pr.dump_stats("cpu_%d.prof"%rank)

if __name__=="__main__":
  run_step5_batch()

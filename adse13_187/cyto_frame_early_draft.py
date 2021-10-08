from __future__ import division, print_function
from time import time
import os, traceback
start_elapse = time()

"""Modified version of the cyto_sim.py script with the following added features:
 - command line parameters following the PHIL + options convention
 - timing statements for weather2.py log scrapping
 - framework for MPI broadcast of big data if and when needed
 - runs on cctbx_project master branch (pickled spectra, for now)
 - Brewster style rank logger
 - Brewster style rank profiles
 - framework to use the exascale API if and when needed
"""

import libtbx.load_env # possibly implicit
from omptbx import omp_get_num_procs
import os,sys

def mcmc_runner_parse_input():
  from LS49.adse13_187.adse13_221.mcmc_runner import generate_phil_scope
  phil_scope = generate_phil_scope()
  params = phil_scope.extract()
  return params

def tsDEPt_one(idx, frame_params):
    print("in tst_run")
    from LS49.adse13_187.adse13_221.mcmc_runner import run as mcmc_run
    params = mcmc_runner_parse_input()
    BERNINA = os.environ.get("BERNINA")
    print("Got",BERNINA)
    params.trusted_mask=os.path.join(BERNINA, "trusted_Py3.mask")
    params.cryst=os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx)
    params.expt = params.cryst
    params.refl=os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx)
    params.output.label="mcmc3"
    params.output.index=idx
    params.model.mosaic_spread.value=0.01
    params.model.Nabc.value=(50,50,15)
    params.model.Nabc.hyperparameter=0.8
    params.model.rot.refine=True
    params.model.cell.covariance=os.path.join(
      BERNINA,"..", "covariance_cytochrome_form.pickle")
    params.simplex.cycles=200
    params.mcmc.cycles=2000
    params.model.plot=False
    print(params.trusted_mask)
    print(params.cryst)
    print(params.expt)
    print(params.refl)
    print(params.model.cell.covariance)
    try:
      mcmc_run(params)
    except Exception as e:
      print("CATCH EXCEPTION",e)

ds1_params="""
ucell_edge_perc=15
ucell_ang_abs=1
method="L-BFGS-B"
spectrum_from_imageset = True
downsamp_spec {
  delta_en = 0.25
}
roi {
  fit_tilt=True
  fit_tilt_using_weights = False
  hotpixel_mask = /global/cscratch1/sd/nksauter/adse13_187/bleededge/work/hopper_help_files/newmask_withbad.pkl
  reject_edge_reflections = False
  pad_shoebox_for_background_estimation=10
}
refiner {
  adu_per_photon = 9.481
  sigma_r=10
}
simulator {
  oversample=4
  structure_factors.mtz_name = /global/cscratch1/sd/nksauter/adse13_187/bleededge/work/hopper_help_files/100shuff.mtz
  structure_factors.mtz_column = "F(+),F(-)"
  beam.size_mm = 0.001
  detector.force_zero_thickness = True
}
init {
  Nabc = 50.000000 50.000000 37.500000
  G = 10.000000
}
mins {
  detz_shift=-1.5
  RotXYZ=[-15,-15,-15]
}
maxs {
  detz_shift = 1.5
  Nabc = 1600 1600 1600
  RotXYZ = 15 15 15
  G = 100000
}
sigmas {
  RotXYZ=[1e-3,1e-3,1e-3]
}
use_restraints = True
betas {
  detz_shift = 1e-08
  ucell = 0.001 0.001
  RotXYZ = 1e-05
  Nabc = 50.000000 50.000000 50.000000
  G = 10.000000
}
centers {
  ucell = 78.61 265.12
  Nabc = 50.000000 50.000000 37.500000
  G = 10.000000
}
fix.detz_shift=True
outdir="."
logging.other_ranks_level="high"
"""


def tst_ds1(idx, frame_params):
    print("in tst_diffbragg_stage_1")

    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz
    from simtbx.diffBragg.phil import hopper_phil
    philz = hopper_phil + philz
    phil_scope_db = parse(philz)
    user_phil = parse(ds1_params)
    working_phil = phil_scope_db.fetch(sources=[user_phil])
    params = working_phil.extract()
    #import logging
    #logging.disable(level=logging.CRITICAL)
    #logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    from simtbx.diffBragg import hopper_utils
    BERNINA = os.environ.get("BERNINA")
    #params.fix.detz_shift = True
    print ("hopper",flush=True)
    try:
      exp, ref, data_modeler, x = hopper_utils.refine(
      exp=os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx),
      ref = os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx),
      params=params, return_modeler=True)
      print("1",flush=True)
      scale, rotX, rotY, rotZ, Na, Nb, Nc,diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c, a,b,c,al,be,ga, detz = hopper_utils.get_param_from_x(x,data_modeler.SIM)
      print("The parameters are\n",scale, rotX, rotY, rotZ, Na, Nb, Nc,diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c, a,b,c,al,be,ga, detz)
    except Exception as e:
      print("CATCH EXCEPTION",e)
      traceback.print_tb(e.__traceback__)

def thin_ds1(idx, frame_params):

    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz
    from simtbx.diffBragg.phil import hopper_phil
    philz = hopper_phil + philz
    phil_scope_db = parse(philz)
    user_phil = parse(ds1_params)
    working_phil = phil_scope_db.fetch(sources=[user_phil])
    params = working_phil.extract()
    from simtbx.diffBragg import hopper_utils
    BERNINA = os.environ.get("BERNINA")

    try:
      exp, ref, data_modeler, x = hopper_utils.refine(
      exp=os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx),
      ref = os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx),
      params=params, return_modeler=True)
      (scale, rotX, rotY, rotZ, Na, Nb, Nc,
       diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c,
       a,b,c,al,be,ga, detz) = hopper_utils.get_param_from_x(x,data_modeler.SIM)
      from dxtbx.model.experiment_list import ExperimentList
      L = ExperimentList()
      L.append(exp)
      L.as_file("ds1_%04d.expt"%idx)
      print("DS1 parameters index %d Na Nb Nc %.3f %.3f %.3f"%(idx,Na,Nb,Nc))

      """
      alldir = os.path.join(os.environ.get("WORK"))
      exptdir = os.path.join(alldir,"2277050")
      allfile = os.path.join(alldir,"allds1.txt")
      with open (allfile,"r") as F:
        lines = F.readlines()
      for line in lines:
        tokens = line.strip().split()
        if int(tokens[3]) == idx:
          Na = float(tokens[7]); Nb = float(tokens[8]); Nc = float(tokens[9]);
          break
      print ("FOUND",idx,end=" ")
      print (Na,Nb,Nc)
      """
      from LS49.adse13_187.adse13_221.basic_runner import run as basic_run
      anaparams = mcmc_runner_parse_input()
      anaparams.trusted_mask=os.path.join(BERNINA, "trusted_Py3.mask")
      anaparams.cryst="ds1_%04d.expt"%idx
      anaparams.expt = os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx)
      anaparams.refl=os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx)
      anaparams.output.enable=False
      anaparams.output.label="mcmc3"
      anaparams.output.index=idx
      anaparams.model.mosaic_spread.value=0.001
      anaparams.model.Nabc.value=(Na,Nb,Nc)
      anaparams.model.Nabc.hyperparameter=0.8
      anaparams.model.rot.refine=True
      anaparams.model.cell.covariance=os.path.join(
      BERNINA,"..", "covariance_cytochrome_form.pickle")
      anaparams.simplex.cycles=200
      anaparams.mcmc.cycles=2000
      anaparams.model.plot=False
      print(anaparams.trusted_mask)
      print(anaparams.cryst)
      print(anaparams.expt)
      print(anaparams.refl)
      print(anaparams.model.cell.covariance)
      basic_run(anaparams)

    except Exception as e:
      print("CATCH EXCEPTION",e)
      traceback.print_exc()


def run_single_job(test_without_mpi=True):
  from LS49.adse13_187.cyto_batch import parse_input as cyto_batch_parse_input
  params,options = cyto_batch_parse_input()

  import datetime
  start_comp = time()

  parcels = list(range(1725,params.N_total))

  print(time(), "delegate parcels")
  os.environ["CCTBX_RECOMMEND_DEVICE"] = "0"

  #os.system("nvidia-smi")
  # client process (requests all the work)
  print("parcels", parcels)
  while len(parcels) > 0:
      idx = parcels[0]
      parcels.remove(idx)
    # server process (does all the work)
      cache_time = time()
      print("idx------start-------->",idx,time(),flush=True)
      thin_ds1(idx,frame_params=params)
      print("idx------finis-------->",idx,
            time(),"elapsed %.3fs"%(time()-cache_time),flush=True)

  print("Overall","at",datetime.datetime.now(),
        "seconds elapsed after srun startup %.3f"%(time()-start_elapse))
  print("Overall","at",datetime.datetime.now(),
        "seconds elapsed after Python imports %.3f"%(time()-start_comp))

def run_batch_job(test_without_mpi=False):
  from LS49.adse13_187.cyto_batch import parse_input as cyto_batch_parse_input
  params,options = cyto_batch_parse_input()
  if params.log.by_rank:
    import io, sys
  if params.log.rank_profile:
    import cProfile
    pr = cProfile.Profile()
    pr.enable()

  if test_without_mpi or params.test_without_mpi:
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
  N_stride = size # total number of worker tasks
  print("hello from rank %d of %d"%(rank,size),"with omp_threads=",omp_get_num_procs())
  import datetime
  start_comp = time()

  if params.log.by_rank:
    expand_dir = os.path.expandvars(params.log.outdir)
    os.makedirs(expand_dir, exist_ok=True)
    log_path = os.path.join(expand_dir,"rank_%d.log"%rank)
    error_path = os.path.join(expand_dir,"rank_%d.err"%rank)
    #print("Rank %d redirecting stdout/stderr to"%rank, log_path, error_path)
    sys.stdout = io.TextIOWrapper(open(log_path,'ab', 0), write_through=True)
    sys.stderr = io.TextIOWrapper(open(error_path,'ab', 0), write_through=True)

  print(rank, time(), "finished with the rank logger, now delgate parcels")
  os.environ["CCTBX_RECOMMEND_DEVICE"] = "%d"%(rank % int(os.environ.get("CCTBX_DEVICE_PER_NODE",1)))
  print("rank", rank, "device", os.environ["CCTBX_RECOMMEND_DEVICE"])
  N_start = int(os.environ.get("N_START",0))

  comm.barrier()
  if rank == 0:
    os.system("nvidia-smi")
    # client process (requests all the work)
    import random
    parcels = list(range(N_start,N_start + params.N_total))
    while len(parcels) > 0:
      idx = parcels[0]    # random.choice(parcels)
      rankreq = comm.recv(source = MPI.ANY_SOURCE)
      print("Sending parcel",idx,"to rank",rankreq)
      comm.send(idx,dest = rankreq)
      parcels.remove(idx)
    # finally send a stop command to each process
    for rankreq in range(size-1):
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send('endrun',dest=rankreq)
  else:
    # server process (does all the work)
    while True:
      # inform the client this worker is ready for an event
      comm.send(rank,dest=0)
      idx = comm.recv(source=0)
      if idx == 'endrun':
        break
      cache_time = time()
      print("idx------start-------->",idx,"rank",rank,time())
      thin_ds1(idx,frame_params=params)
      print("idx------finis-------->",idx,
            "rank",rank,time(),"elapsed %.3fs"%(time()-cache_time))
  comm.barrier()

  print("Overall rank",rank,"at",datetime.datetime.now(),
        "seconds elapsed after srun startup %.3f"%(time()-start_elapse))
  print("Overall rank",rank,"at",datetime.datetime.now(),
        "seconds elapsed after Python imports %.3f"%(time()-start_comp))
  if params.log.rank_profile:
    pr.disable()
    pr.dump_stats("cpu_%d.prof"%rank)

if __name__=="__main__":
  run_batch_job()
  #run_single_job()

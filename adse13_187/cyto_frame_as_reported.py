from __future__ import division, print_function
from time import time
import os, traceback
import numpy as np
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

def thin_ds1(idx, frame_params):

    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz
    from simtbx.diffBragg.phil import hopper_phil
    philz = hopper_phil + philz
    phil_scope_db = parse(philz)
    from LS49.adse13_187.report_versions import ds1_params_v4 as ds1_params
    user_phil = parse(ds1_params)
    working_phil = phil_scope_db.fetch(sources=[user_phil])
    diffbragg_params = working_phil.extract()
    from simtbx.diffBragg import hopper_utils
    BERNINA = os.environ.get("BERNINA")
    anaparams = mcmc_runner_parse_input()
    if frame_params.use_diffuse_models:
      diffbragg_params.use_diffuse_models = True
      diffbragg_params.method = "Nelder-Mead"
      diffbragg_params.init.diffuse_sigma = 1,1,1
      diffbragg_params.init.diffuse_gamma = 100,100,100
      diffbragg_params.betas.G = 1000.
      diffbragg_params.fix.diffuse_gamma = False
      diffbragg_params.fix.diffuse_sigma = False

    try:
      exp, ref, data_modeler, x = hopper_utils.refine(
      exp=os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx),
      ref = os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx),
      params=diffbragg_params, return_modeler=True)

      "This subsection writes out the diffBragg-calculated model image"
      """
      from simtbx.command_line.hopper import save_to_pandas
      stage1_df = save_to_pandas(x, data_modeler.SIM, os.path.join(BERNINA, "split_cs", "split_%04d.expt"%idx),
              diffbragg_params,data_modeler.E, 0, os.path.join(BERNINA, "split2b" , "split_%04d.refl"%idx), None)
      stage1_df.to_pickle("stage1_df_%04d.pickle"%idx)
      np.save("modeler_%04d"%idx, data_modeler)
      hopper_utils.write_SIM_logs(data_modeler.SIM, log="state_fname_%04d.txt"%idx, lam="spectra_fname_%04d.txt"%idx)
      """

      (scale, rotX, rotY, rotZ, Na, Nb, Nc,
       diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c,
       a,b,c,al,be,ga, detz) = hopper_utils.get_param_from_x(x,data_modeler.SIM)
      from dxtbx.model.experiment_list import ExperimentList
      L = ExperimentList()
      L.append(exp)
      L.as_file("ds1_%04d.expt"%idx)
      print("DS1 parameters index %d Na Nb Nc %.3f %.3f %.3f"%(idx,Na,Nb,Nc))

      from LS49.adse13_187.adse13_221.basic_runner import run as basic_run
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
      if frame_params.use_diffuse_models:
        basic_run(anaparams,use_diffuse=(diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c))
      else:
        basic_run(anaparams)

    except Exception as e:
      print("CATCH EXCEPTION",e)
      traceback.print_exc()

def run_single_job(test_without_mpi=True):
  from LS49.adse13_187.cyto_batch import parse_input as cyto_batch_parse_input
  params,options = cyto_batch_parse_input()

  import datetime
  start_comp = time()

  parcels = list(range(int(os.environ.get("N_START"),0),params.N_total))

  print(time(), "delegate parcels")
  os.environ["CCTBX_RECOMMEND_DEVICE"] = "0"

  #os.system("nvidia-smi")
  # client process (requests all the work)
  print("parcels", parcels)

  # with this change, I would hope you no longer need to edit out the MAIN_LOGGER lines
  import logging
  dblogger = logging.getLogger("diffBragg.main")
  handler = logging.StreamHandler()
  dblogger.setLevel(logging.DEBUG)
  handler.setLevel(logging.DEBUG)
  # to write less output, increase the logging level, these are just pointers to numbers,
  # logging.INFO will write less output, logging.CRITICAL will mostly silience diffBragg
  dblogger.addHandler(handler) # In the hopper_utils library, logging will use this stream

  while len(parcels) > 0:
      idx = parcels[0]
      parcels.remove(idx)
    # server process (does all the work)
      cache_time = time()
      print("idx------start-------->",idx,time(),flush=True)
      #from LS49.adse13_187.cyto_frame_early_draft import tsDEPt_one
      #tsDEPt_one(idx,frame_params=params)
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

    # ON LOGGING,
    # with this change, I would hope you no longer need to edit out the MAIN_LOGGER lines
    import logging
    dblogger = logging.getLogger("diffBragg.main")
    # we can change 'main' to e.g. 'diffBragg.main', but I would
    # need to edit hopper_utils.py to reflect that, just a one-line fix
    handler = logging.StreamHandler(stream=sys.stdout)
    # handler.setLevel(logging.DEBUG) # not needed
    dblogger.setLevel(logging.DEBUG)
    # to write less output, increase the logging level, these are just pointers to numbers,
    # logging.INFO will write less output, logging.CRITICAL will mostly silience diffBragg
    dblogger.addHandler(handler) # Now, in the hopper_utils library, logging shoud use this stream

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
      #from LS49.adse13_187.cyto_frame_early_draft import tsDEPt_one
      #tsDEPt_one(idx,frame_params=params)
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

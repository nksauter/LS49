from __future__ import division, print_function
from time import time
start_elapse = time()

import libtbx.load_env # possibly implicit
from omptbx import omp_get_num_procs
import os,sys

def parse_input():
  from iotbx.phil import parse
  master_phil="""
    N_total = 75
      .type = int
      .help = number of events to simulate
    log {
      outdir = .
        .type = path
        .help = Use "/mnt/bb/${USER}" for Summit NVME burst buffer
      by_rank = True
        .type = bool
        .help = Brewster-style split of logs for each rank
      rank_profile = False
        .type = bool
        .help = create a cProfile output for each rank
    }
    devices_per_node = 1
      .type = int
      .help = always 1 per Summit resource group, either 1 or 8 for Cori GPU
    use_exascale_api = True
      .type = bool
      .help = aim for 3 second image turnaround
  """
  phil_scope = parse(master_phil)
  # The script usage
  import libtbx.load_env # implicit import
  help_message = '''Jungfrau/cytochrome simuation.'''
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

def tst_one(i_exp,spectra,Fmerge,gpu_channels_singleton,rank,params):
    from simtbx.nanoBragg import utils
    from dxtbx.model.experiment_list import ExperimentListFactory
    import numpy as np

    print("Experiment %d" % i_exp, flush=True)
    sys.stdout.flush()

    save_data_too = True
    outfile = "boop_%d.hdf5" % i_exp
    experiment_file = "/global/cfs/cdirs/m3562/der/run795/top_%d.expt" % i_exp
    refl_file = "/global/cfs/cdirs/m3562/der/run795/top_%d.refl" % i_exp
    cuda = True  # False  # whether to use cuda
    omp = False
    ngpu_on_node = 1 # 8  # number of available GPUs
    mosaic_spread = 0.07  # degrees
    mosaic_spread_samples = 500 # 30  # 50  # number of mosaic blocks sampling mosaicity
    Ncells_abc = 30, 30, 10  # medians from best stage1
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    plot_spec = False  # plot the downsample spectra before simulating
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    panel_list = None  # integer list of panels, usefule for debugging
    rois_only = False  # only set True if you are running openMP, or CPU-only (i.e. not for GPU)
    include_background = True  # default is to add water background 100 mm thick
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    # <><><><><><><><><><><><><><><><>

    El = ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=save_data_too)
    exper = El[0]
    crystal = exper.crystal
    detector = exper.detector
    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)

    beam = exper.beam

    if flat:
        assert detector[0].get_thickness() == 0

    if panel_list is None:
        panel_list = list(range(len(detector)))

    energies, weights = spectra[i_exp]

    pids_for_rank = panel_list
    device_Id = gpu_channels_singleton.get_deviceID()
    print("Rank %d will use device %d" % (rank, device_Id))
    show_params = (rank == 0)  # False
    time_panels = (rank == 0)

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = utils.ENERGY_CONV / mn_energy

    if params.use_exascale_api:
      utils.multipanel_sim(
      CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
      energies=energies, fluxes=weights, Famp=Fmerge,
      cuda=cuda, device_Id=device_Id,
      oversample=oversample, Ncells_abc=Ncells_abc, verbose=verbose,
      time_panels=time_panels, show_params=show_params, spot_scale_override=spot_scale,
      mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread, beamsize_mm=beamsize_mm,
      )
      exit("23")

      #optional background
      backgrounds = {pid: None for pid in panel_list}
      if include_background:
        backgrounds = {pid: utils.sim_background( # default is for water
                detector, beam, wavelengths=[mn_wave], wavelength_weights=[1],
                total_flux=total_flux,
                pidx=pid, beam_size_mm=beamsize_mm, sample_thick_mm=0.5) # 0.1)
            for pid in pids_for_rank}

      include_noise=False

    #optional background
    backgrounds = {pid: None for pid in panel_list}
    if include_background:
        backgrounds = {pid: utils.sim_background( # default is for water
                detector, beam, wavelengths=[mn_wave], wavelength_weights=[1],
                total_flux=total_flux,
                pidx=pid, beam_size_mm=beamsize_mm, sample_thick_mm=0.5) # 0.1)
            for pid in pids_for_rank}

    pid_and_pdata = utils.flexBeam_sim_colors(
      CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
      energies=energies, fluxes=weights, Famp=Fmerge,
      pids=pids_for_rank, cuda=cuda, device_Id=device_Id,
      oversample=oversample, Ncells_abc=Ncells_abc, verbose=verbose,
      time_panels=time_panels, show_params=show_params, spot_scale_override=spot_scale,
      mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread, beamsize_mm=beamsize_mm,
      background_raw_pixels=backgrounds, include_noise=False, rois_perpanel=None)

    pid_and_pdata = sorted(pid_and_pdata, key=lambda x: x[0])
    _, pdata = zip(*pid_and_pdata)

    # pdata is a list of 256 2D numpy arrays, now.
    exit("44")

    if len(panel_list) != len(detector):
        print("Cant save partial detector image, exiting..")
        exit()
        #from dxtbx.model import Detector
        #new_det = Detector()
        #for pid in panel_list:
        #    new_det.add_panel(detector[pid])
        #detector = new_det
    if save_data_too:
        data = exper.imageset.get_raw_data(0)

    tsave = time()
    pdata = np.array(pdata) # now pdata is a numpy array of shape 256,254,254
    img_sh = pdata.shape
    num_output_images = 1 + int(save_data_too)
    print("Saving output data of shape", img_sh)
    beam_dict = beam.to_dict()
    det_dict = detector.to_dict()
    try:
      beam_dict.pop("spectrum_energies")
      beam_dict.pop("spectrum_weights")
    except Exception: pass
    with utils.H5AttributeGeomWriter(outfile, image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict,
                                detector_and_beam_are_dicts=True) as writer:
        writer.add_image(pdata)

        if save_data_too:
            data = [data[pid].as_numpy_array() for pid in panel_list]
            writer.add_image(data)

    tsave = time() - tsave
    print("Saved output to file %s. Saving took %.4f sec" % (outfile, tsave, ))


def run_batch_job(test_without_mpi=False):
  params,options = parse_input()
  if params.log.by_rank:
    import io, sys
  if params.log.rank_profile:
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
  N_stride = size # total number of worker tasks
  print("hello from rank %d of %d"%(rank,size),"with omp_threads=",omp_get_num_procs())
  import datetime
  start_comp = time()

  print(rank, time(),
    "finished with the calculation of channels, now construct single broadcast")

  if rank == 0:
    print("Rank 0 time", datetime.datetime.now())

    spectrum_dict = {}
    with open("../test.pickle","rb") as F:
      for i_exp in range(75):
        import pickle
        i_exp_p, energies, weights = pickle.load(F)
        assert i_exp == i_exp_p
        spectrum_dict[i_exp] = (energies, weights)

    from iotbx.reflection_file_reader import any_reflection_file
    merge_file = "/global/cfs/cdirs/m3562/der/cyto_init_merge.mtz"
    Fmerge = any_reflection_file(merge_file).as_miller_arrays()[0].as_amplitude_array()
    if comm.rank == 0:
        print("Fmerge min/max = %f / %f" % (min(Fmerge.data()), max(Fmerge.data())))

    transmitted_info = dict(spectra = spectrum_dict,
                            amplitudes = Fmerge,
                            )
  else:
    transmitted_info = None
  transmitted_info = comm.bcast(transmitted_info, root = 0)
  comm.barrier()
  parcels = list(range(rank,params.N_total,N_stride))

  print(rank, time(), "finished with single broadcast, now set up the rank logger")

  if params.log.by_rank:
    expand_dir = os.path.expandvars(params.log.outdir)
    log_path = os.path.join(expand_dir,"rank_%d.log"%rank)
    error_path = os.path.join(expand_dir,"rank_%d.err"%rank)
    #print("Rank %d redirecting stdout/stderr to"%rank, log_path, error_path)
    sys.stdout = io.TextIOWrapper(open(log_path,'ab', 0), write_through=True)
    sys.stderr = io.TextIOWrapper(open(error_path,'ab', 0), write_through=True)

  print(rank, time(), "finished with the rank logger, now construct the GPU cache container")

  from simtbx.nanoBragg import gpu_energy_channels
  gpu_channels_singleton = gpu_energy_channels (
    deviceId = rank % params.devices_per_node )
    # singleton will instantiate, regardless of cuda, device count, or exascale API

  comm.barrier()
  import random
  while len(parcels)>0:
    idx = random.choice(parcels)
    cache_time = time()
    print("idx------start-------->",idx,"rank",rank,time())
    # if rank==0: os.system("nvidia-smi")
    tst_one(i_exp=idx,spectra=transmitted_info["spectra"],
        Fmerge=transmitted_info["amplitudes"],
        gpu_channels_singleton=gpu_channels_singleton,
        rank=rank,params=params
    )
    parcels.remove(idx)
    print("idx------finis-------->",idx,"rank",rank,time(),"elapsed",time()-cache_time)
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

from __future__ import division, print_function
from time import time
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
from xfel.merging.application.utils.memory_usage import get_memory_usage
import os,sys

from simtbx.nanoBragg.tst_gauss_argchk import water
# water is a flex.vec2_double, with structure factor background vs. sin-theta over lambda

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
        .help = Set to False for IPython debugging
      rank_profile = False
        .type = bool
        .help = create a cProfile output for each rank
    }
    test_without_mpi = False
      .type = bool
      .help = True forces the one-rank MPI emulator
    devices_per_node = 1
      .type = int
      .help = always 1 per Summit resource group, either 1 or 8 for Cori GPU
    use_exascale_api = True
      .type = bool
      .help = aim for 3 second image turnaround
    test_pixel_congruency = False
      .type = bool
      .help = ensure the new /old style agree at per-pixel level
    include_background = True
      .type = bool
      .help = whether to add background to model
    write_output = False
      .type = bool
      .help = whether to write an output image file (hdf5)
    write_experimental_data = False
      .type = bool
      .help = if hdf5 file is written, also include a frame giving the experimental data
      .help = needs to be a command line argument in case the experimental data are unavailable
    mosaic_spread_samples = 500
      .type = int
      .help = granularity of mosaic rotation, double it to find number of umats
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

def multipanel_sim(
  CRYSTAL, DETECTOR, BEAM, Famp, energies, fluxes,
  background_wavelengths, background_wavelength_weights,
  background_total_flux, background_sample_thick_mm,
  density_gcm3=1, molecular_weight=18,
  cuda=False, oversample=0, Ncells_abc=(50, 50, 50),
  mos_dom=1, mos_spread=0, beamsize_mm=0.001,
  crystal_size_mm=0.01,
  verbose=0, default_F=0, interpolate=0, profile="gauss",
  spot_scale_override=None, show_params=False, time_panels=False,
  add_water = False, add_air=False, water_path_mm=0.005, air_path_mm=0,
  adc_offset=0, readout_noise=3, psf_fwhm=0, gain=1, mosaicity_random_seeds=None, include_background=True):

  from simtbx.nanoBragg.nanoBragg_beam import NBbeam
  from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
  from simtbx.nanoBragg.sim_data import SimData
  from simtbx.nanoBragg.utils import get_xray_beams
  from scitbx.array_family import flex
  from scipy import constants
  import numpy as np
  ENERGY_CONV = 10000000000.0 * constants.c * constants.h / constants.electron_volt
  assert cuda # disallow the default

  nbBeam = NBbeam()
  nbBeam.size_mm = beamsize_mm
  nbBeam.unit_s0 = BEAM.get_unit_s0()
  wavelengths = ENERGY_CONV / np.array(energies)
  nbBeam.spectrum = list(zip(wavelengths, fluxes))

  nbCrystal = NBcrystal()
  nbCrystal.dxtbx_crystal = CRYSTAL
  #nbCrystal.miller_array = None # use the gpu_channels_singleton mechanism instead
  nbCrystal.Ncells_abc = Ncells_abc
  nbCrystal.symbol = CRYSTAL.get_space_group().info().type().lookup_symbol()
  nbCrystal.thick_mm = crystal_size_mm
  nbCrystal.xtal_shape = profile
  nbCrystal.n_mos_domains = mos_dom
  nbCrystal.mos_spread_deg = mos_spread

  pid = 0 # remove the loop, use C++ iteration over detector panels
  use_exascale_api = True
  if use_exascale_api:
    tinit = time()
    S = SimData()
    S.detector = DETECTOR
    S.beam = nbBeam
    S.crystal = nbCrystal
    S.panel_id = pid
    S.add_air = add_air
    S.air_path_mm = air_path_mm
    S.add_water = add_water
    S.water_path_mm = water_path_mm
    S.readout_noise = readout_noise
    S.gain = gain
    S.psf_fwhm = psf_fwhm
    S.include_noise = False

    if mosaicity_random_seeds is not None:
      S.mosaic_seeds = mosaicity_random_seeds

    S.instantiate_nanoBragg(verbose=verbose, oversample=oversample, interpolate=interpolate,
      device_Id=Famp.get_deviceID(),default_F=default_F, adc_offset=adc_offset)

    SIM = S.D # the nanoBragg instance
    assert Famp.get_deviceID()==SIM.device_Id
    if spot_scale_override is not None:
      SIM.spot_scale = spot_scale_override
    assert Famp.get_nchannels() == 1 # non-anomalous scenario

    from simtbx.gpu import exascale_api
    gpu_simulation = exascale_api(nanoBragg = SIM)
    gpu_simulation.allocate_cuda() # presumably done once for each image

    from simtbx.gpu import gpu_detector as gpud
    gpu_detector = gpud(deviceId=SIM.device_Id, detector=DETECTOR,
                        beam=BEAM)
    gpu_detector.each_image_allocate_cuda()

    # revisit the allocate cuda for overlap with detector, sync up please
    x = 0 # only one energy channel
    #SIM.flux = background_total_flux
      #SIM.flux = self.flux[x]
      #SIM.wavelength_A = self.wavlen[x]
    from libtbx.development.timers import Profiler
    P = Profiler("from gpu amplitudes cuda")
    gpu_simulation.add_energy_channel_from_gpu_amplitudes_cuda(
      x, Famp, gpu_detector)
    TIME_BRAGG = time()-P.start_el
    del P
#now in position to implement the detector panel loop

    per_image_scale_factor = 1./len(energies)
    gpu_detector.scale_in_place_cuda(per_image_scale_factor) # apply scale directly on GPU

    if include_background:
      t_bkgrd_start = time()
      SIM.beamsize_mm = beamsize_mm

      wavelength_weights = np.array(background_wavelength_weights)
      weights = wavelength_weights / wavelength_weights.sum() * background_total_flux
      spectrum = list(zip(background_wavelengths, weights))
      xray_beams = get_xray_beams(spectrum, BEAM)
      SIM.xray_beams = xray_beams
      SIM.Fbg_vs_stol = water
      SIM.flux=sum(weights)
      SIM.amorphous_sample_thick_mm = background_sample_thick_mm
      SIM.amorphous_density_gcm3 = density_gcm3
      SIM.amorphous_molecular_weight_Da = molecular_weight
      gpu_simulation.add_background_cuda(gpu_detector)
      TIME_BG = time()-t_bkgrd_start

    packed_numpy = gpu_detector.get_raw_pixels_cuda().as_numpy_array()
    gpu_detector.each_image_free_cuda()

    return packed_numpy, TIME_BG, TIME_BRAGG

def tst_one(i_exp,spectra,Fmerge,gpu_channels_singleton,rank,params):
    from simtbx.nanoBragg import utils
    from dxtbx.model.experiment_list import ExperimentListFactory
    import numpy as np

    print("Experiment %d" % i_exp, flush=True)
    sys.stdout.flush()

    outfile = "boop_%d.hdf5" % i_exp
    from LS49.adse13_187.case_data import retrieve_from_repo
    experiment_file = retrieve_from_repo(i_exp)
    # Not used # refl_file = "/global/cfs/cdirs/m3562/der/run795/top_%d.refl" % i_exp
    cuda = True  # False  # whether to use cuda
    omp = False
    ngpu_on_node = 1 # 8  # number of available GPUs
    mosaic_spread = 0.07  # degrees
    mosaic_spread_samples = params.mosaic_spread_samples # number of mosaic blocks sampling mosaicity
    Ncells_abc = 30, 30, 10  # medians from best stage1
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    plot_spec = False  # plot the downsample spectra before simulating
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    panel_list = None  # integer list of panels, usefule for debugging
    rois_only = False  # only set True if you are running openMP, or CPU-only (i.e. not for GPU)
    include_background = params.include_background   # default is to add water background 100 mm thick
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    #<><><><><><><><>
    # XXX new code
    El = ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=True)
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

    # XXX new code
    spec = exper.imageset.get_spectrum(0)
    energies_raw, weights_raw = spec.get_energies_eV().as_numpy_array(), \
                                spec.get_weights().as_numpy_array()
    energies, weights = utils.downsample_spectrum(energies_raw, weights_raw, method=1, total_flux=total_flux,
                                                  ev_width=ev_res)

    if flat:
        assert detector[0].get_thickness() == 0

    if panel_list is None:
        panel_list = list(range(len(detector)))

    pids_for_rank = panel_list
    device_Id = 0
    if gpu_channels_singleton is not None:
      device_Id = gpu_channels_singleton.get_deviceID()

    print("Rank %d will use device %d" % (rank, device_Id))
    show_params = False
    time_panels = (rank == 0)

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = utils.ENERGY_CONV / mn_energy

    if params.use_exascale_api:
      BEG=time()
      print (gpu_channels_singleton.get_deviceID(),"device")
      Famp_is_uninitialized = ( gpu_channels_singleton.get_nchannels() == 0 ) # uninitialized
      if Famp_is_uninitialized:
        F_P1 = Fmerge.expand_to_p1()
        for x in range(1):  # in this scenario, amplitudes are independent of lambda
          gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
      assert gpu_channels_singleton.get_nchannels() == 1

      JF16M_numpy_array, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
        Famp = gpu_channels_singleton,
        energies=list(energies), fluxes=list(weights),
        background_wavelengths=[mn_wave], background_wavelength_weights=[1],
        background_total_flux=total_flux,background_sample_thick_mm=0.5,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
        beamsize_mm=beamsize_mm,show_params=show_params,
        time_panels=time_panels, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=include_background)
      TIME_EXA = time()-BEG
      print ("Exascale time",TIME_EXA)
      if params.write_experimental_data:
        data = exper.imageset.get_raw_data(0)

      tsave = time()
      img_sh = JF16M_numpy_array.shape
      assert img_sh == (256,254,254)
      num_output_images = 1 + int(params.write_experimental_data)
      print("Saving exascale output data of shape", img_sh)
      beam_dict = beam.to_dict()
      det_dict = detector.to_dict()
      try:
        beam_dict.pop("spectrum_energies")
        beam_dict.pop("spectrum_weights")
      except Exception: pass
# XXX no longer have two separate files
      if params.write_output:
       with utils.H5AttributeGeomWriter("exap_%d.hdf5"%i_exp,
                                image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict,
                                detector_and_beam_are_dicts=True) as writer:
        writer.add_image(JF16M_numpy_array)

        if params.write_experimental_data:
            data = [data[pid].as_numpy_array() for pid in panel_list]
            writer.add_image(data)

       tsave = time() - tsave
       print("Saved output to file %s. Saving took %.4f sec" % ("exap_%d.hdf5"%i_exp, tsave, ))

    BEG2 = time()
    #optional background
    TIME_BG2 = time()
    backgrounds = {pid: None for pid in panel_list}
    if include_background:
        backgrounds = {pid: utils.sim_background( # default is for water
                detector, beam, wavelengths=[mn_wave], wavelength_weights=[1],
                total_flux=total_flux, Fbg_vs_stol=water,
                pidx=pid, beam_size_mm=beamsize_mm, sample_thick_mm=0.5)
            for pid in pids_for_rank}
    TIME_BG2 = time()-TIME_BG2

    TIME_BRAGG2 = time()
    pid_and_pdata = utils.flexBeam_sim_colors(
      CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
      energies=list(energies), fluxes=list(weights), Famp=Fmerge,
      pids=pids_for_rank, cuda=cuda, device_Id=device_Id,
      oversample=oversample, Ncells_abc=Ncells_abc, verbose=verbose,
      time_panels=time_panels, show_params=show_params, spot_scale_override=spot_scale,
      mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread, beamsize_mm=beamsize_mm,
      background_raw_pixels=backgrounds, include_noise=False, rois_perpanel=None)
    TIME_BRAGG2 = time()-TIME_BRAGG2
    pid_and_pdata = sorted(pid_and_pdata, key=lambda x: x[0])
    _, pdata = zip(*pid_and_pdata)
    TIME_VINTAGE = time()-BEG2

    print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("\tBreakdown:")
    if params.use_exascale_api:
        print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    print("\t\tVintage:  time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG2, TIME_BRAGG2, TIME_VINTAGE))
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")

    if params.test_pixel_congruency and params.use_exascale_api:
      abs_diff = np.abs(np.array(pdata) - JF16M_numpy_array).max()
      assert np.allclose(pdata, JF16M_numpy_array), "max per-pixel difference: %f photons"%abs_diff
      print("pixel congruency: OK!")

    # pdata is a list of 256 2D numpy arrays, now.

    if len(panel_list) != len(detector):
        print("Cant save partial detector image, exiting..")
        exit()
        #from dxtbx.model import Detector
        #new_det = Detector()
        #for pid in panel_list:
        #    new_det.add_panel(detector[pid])
        #detector = new_det
    if params.write_experimental_data:
        data = exper.imageset.get_raw_data(0)

    tsave = time()
    pdata = np.array(pdata) # now pdata is a numpy array of shape 256,254,254
    img_sh = pdata.shape
    num_output_images = 3 + int(params.write_experimental_data)
    print("BOOPZ: Rank=%d ; i_exp=%d, RAM usage=%f" % (rank, i_exp,get_memory_usage()/1e6 ))
    beam_dict = beam.to_dict()
    det_dict = detector.to_dict()
    try:
      beam_dict.pop("spectrum_energies")
      beam_dict.pop("spectrum_weights")
    except Exception: pass
    if params.write_output:
      print("Saving output data of shape", img_sh)
      with utils.H5AttributeGeomWriter(outfile, image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict,
                                detector_and_beam_are_dicts=True) as writer:
        writer.add_image(JF16M_numpy_array/pdata)
        writer.add_image(JF16M_numpy_array)
        writer.add_image(pdata)

        if params.write_experimental_data:
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

  print(rank, time(),
    "finished with the calculation of channels, now construct single broadcast")

  if rank == 0:
    print("Rank 0 time", datetime.datetime.now())

    spectrum_dict = {}

    from iotbx.reflection_file_reader import any_reflection_file
    from LS49 import ls49_big_data
    merge_file = os.path.join(ls49_big_data,"adse13_228","cyto_init_merge.mtz")
    Fmerge = any_reflection_file(merge_file).as_miller_arrays()[0].as_amplitude_array()

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
    os.makedirs(expand_dir, exist_ok=True)
    log_path = os.path.join(expand_dir,"rank_%d.log"%rank)
    error_path = os.path.join(expand_dir,"rank_%d.err"%rank)
    #print("Rank %d redirecting stdout/stderr to"%rank, log_path, error_path)
    sys.stdout = io.TextIOWrapper(open(log_path,'ab', 0), write_through=True)
    sys.stderr = io.TextIOWrapper(open(error_path,'ab', 0), write_through=True)

  print(rank, time(), "finished with the rank logger, now construct the GPU cache container")

  try:
    from simtbx.gpu import gpu_energy_channels
    gpu_channels_singleton = gpu_energy_channels (
      deviceId = rank % params.devices_per_node )
      # singleton will instantiate, regardless of cuda, device count, or exascale API
  except ImportError:
    gpu_channels_singleton = None
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

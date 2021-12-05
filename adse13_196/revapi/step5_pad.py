from __future__ import division,print_function
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal
import math
import scitbx
import os
from libtbx.development.timers import Profiler
from simtbx.gpu import exascale_api, gpu_detector as gpud
from dxtbx.model.detector import DetectorFactory

big_data = "." # directory location for reference files
def full_path(filename):
  return os.path.join(big_data,filename)

def data():
  from LS49.sim.fdp_plot import george_sherrell
  return dict(
    pdb_lines = open(full_path("1m2a.pdb"),"r").read(),
    Fe_oxidized_model = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out")),
    Fe_reduced_model = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out")),
    Fe_metallic_model = george_sherrell(full_path("data_sherrell/Fe_fake.dat"))
  )

def raw_to_pickle(raw_pixels, fileout):
  from six.moves import cPickle as pickle
  with open(fileout, "wb") as F:
    pickle.dump(raw_pixels, F)

from LS49.sim.step4_pad import microcrystal # implicit import

"""Changes in Step4K relative to Step4
k_sol is now hard-coded to 0.435 instead of 0.35 to reduce solvent contrast
Simulation out to 1.7 Angstrom resolution
Use 3K x 3K detector, not 2K x 2K to insure complete coverage at 2.1 A inscribed circle
Put the ROTMAT orientation in the file header, instead of in log file
Put sequence number and MPI rank in header
Write gz compressed file
Avoid computing the image if already on disk

Changes in step5 relative to step4K
change to wavelength 7120 eV
Change pathlength from 1 um to 4 um. Low angle spots may be overloaded, but hope for higher-fidelity at 2.3 Angstrom resolution.
Calculate the scattering factors separately for each wavelength channel in the incident spectrum, using gen_fmodel class.
Identify which iron atoms are reduced and oxidized in "red-fd" 1m2a:
  buried: Cys22-Fe1 A201-Cys9, Fe1 B202, higher energy, oxidized, Fe(III), ferric
 surface: Cys59-Fe2 A201-Cys55, Fe2 B202, lower energy, reduced, Fe(II), ferrous
Reset fp and fdp for identified scatterers
"""

"""Prospective plan
Simulate an anomalous dataset, 7150 eV, but using fixed f',f". 50,000 images. Solve by SAD. (step4K)
*** Plug in the f' & f" as a function of wavelength.  Rerun the simulation at 7120 eV. (step5)
  Sort images into normalized-intensity wavelength bins and calculate maps
*** Put a detector at far distance, and one nearby, process both simultaneously.
  How to modify nanoBragg so it takes an offset Jungfrau, CSPAD, or Rayonix
  Make sure we can index both at the same time.
  Observe the energy positioning of Bragg spots, as function of deltapsi.
*** simulate polarized spectroscopy
"""

"""dials integration
dxtbx.print_header step4_000001.img
dials.import step4_000001.img # checks import
dials.estimate_gain datablock.json
dials.find_spots datablock.json threshold.dispersion.gain=1.47 filter.min_spot_size=2
dials.image_viewer datablock.json strong.pickle
dials.stills_process step4_00000[0-2].img threshold.dispersion.gain=1.47 filter.min_spot_size=2 indexing.known_symmetry.unit_cell=67.200,59.800,47.200,90.00,110.30,90.00 indexing.known_symmetry.space_group=C2 mp.nproc=60
dials.image_viewer idx-step4_000000_integrated_experiments.json idx-step4_000000_integrated.pickle
"""
def write_safe(fname):
  # make sure file or compressed file is not already on disk
  return (not os.path.isfile(fname)) and (not os.path.isfile(fname+".gz"))

add_background_algorithm = str(os.environ.get("ADD_BACKGROUND_ALGORITHM","jh"))
assert add_background_algorithm in ["jh","sort_stable","cuda"]

def basic_detector(): # from development, not actually used here but useful for reference
  # make a detector panel
  # monolithic camera description
  print("Make a dxtbx detector")
  detdist = 141.7
  pixsize = 0.11
  im_shape = 3000, 3000
  det_descr = {'panels':
               [{'fast_axis': (1.0, 0.0, 0.0),
                 'slow_axis': (0.0, -1.0, 0.0),
                 'gain': 1.0,
                 'identifier': '',
                 'image_size': im_shape,
                 'mask': [],
                 'material': '',
                 'mu': 0.0,
                 'name': 'Panel',
                 'origin': (-im_shape[0]*pixsize/2., im_shape[1]*pixsize/2., -detdist),
                 'pedestal': 0.0,
                 'pixel_size': (pixsize, pixsize),
                 'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                 'raw_image_offset': (0, 0),
                 'thickness': 0.0,
                 'trusted_range': (-1e7, 1e7),
                 'type': ''}]}
  return DetectorFactory.from_dict(det_descr)

from LS49.sim.debug_utils import channel_extractor
CHDBG_singleton = channel_extractor()

def run_sim2smv(prefix,crystal,spectra,rotation,rank,gpu_channels_singleton,params,
                quick=False,save_bragg=False,sfall_channels=None):
  smv_fileout = prefix + ".img"
  burst_buffer_expand_dir = os.path.expandvars(params.logger.outdir)
  burst_buffer_fileout = os.path.join(burst_buffer_expand_dir,smv_fileout)
  reference_fileout = os.path.join(".",smv_fileout)
  if not quick:
    if not write_safe(reference_fileout):
      print("File %s already exists, skipping in rank %d"%(reference_fileout,rank))
      return

  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = next(spectra) # list of lambdas, list of fluxes, average wavelength
  assert wavelength_A > 0

  # use crystal structure to initialize Fhkl array
  N = crystal.number_of_cells(sfall_channels[0].unit_cell())

  #SIM = nanoBragg(detpixels_slowfast=(2000,2000),pixel_size_mm=0.11,Ncells_abc=(5,5,5),verbose=0)
  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    # workaround for problem with wavelength array, specify it separately in constructor.
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 0 # Do not offset by 40
  SIM.adc_offset_adu = 10 # Do not offset by 40
  import sys
  if len(sys.argv)>2:
    SIM.seed = -int(sys.argv[2])
  if len(sys.argv)>1:
    if sys.argv[1]=="random" : SIM.randomize_orientation()
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
                               # mosaic_domains setter MUST come after mosaic_spread_deg setter
  SIM.mosaic_domains = int(os.environ.get("MOS_DOM","25"))
  print ("MOSAIC",SIM.mosaic_domains)
  SIM.distance_mm=141.7

  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  scitbx.random.set_random_seed(1234)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=SIM.mosaic_spread_deg * math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(SIM.mosaic_domains)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )
  SIM.set_mosaic_blocks(UMAT_nm)

  #SIM.detector_thick_mm = 0.5 # = 0 for Rayonix
  #SIM.detector_thicksteps = 1 # should default to 1 for Rayonix, but set to 5 for CSPAD
  #SIM.detector_attenuation_length_mm = default is silicon

  # get same noise each time this test is run
  SIM.seed = 1
  SIM.oversample=1
  SIM.wavelength_A = wavelength_A
  SIM.polarization=1
  # this will become F000, marking the beam center
  SIM.default_F=0
  SIM.Fhkl=sfall_channels[0] # instead of sfall_main
  Amatrix_rot = (rotation *
             sqr(sfall_channels[0].unit_cell().orthogonalization_matrix())).transpose()

  SIM.Amatrix_RUB = Amatrix_rot
  #workaround for failing init_cell, use custom written Amatrix setter
  print("unit_cell_Adeg=",SIM.unit_cell_Adeg)
  print("unit_cell_tuple=",SIM.unit_cell_tuple)
  Amat = sqr(SIM.Amatrix).transpose() # recovered Amatrix from SIM
  from cctbx import crystal_orientation
  Ori = crystal_orientation.crystal_orientation(Amat, crystal_orientation.basis_type.reciprocal)

  # fastest option, least realistic
  SIM.xtal_shape=shapetype.Gauss_argchk # both crystal & RLP are Gaussian
  # only really useful for long runs
  SIM.progress_meter=False
  # prints out value of one pixel only.  will not render full image!
  # flux is always in photons/s
  SIM.flux=1e12
  SIM.exposure_s=1.0 # so total fluence is e12
  # assumes round beam
  SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
  temp=SIM.Ncells_abc
  SIM.Ncells_abc=temp

  # rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
  water_bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.18,7.32),(0.2,6.75),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
  assert [a[0] for a in water_bg] == sorted([a[0] for a in water_bg])
  # rough approximation to air
  air_bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),(0.5,4.22)])
  assert [a[0] for a in air_bg] == sorted([a[0] for a in air_bg])

  # simulated crystal is only 125 unit cells (25 nm wide)
  # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
  SIM.raw_pixels *= crystal.domains_per_crystal; # must calculate the correct scale!

  QQ = Profiler("nanoBragg Bragg spots rank %d"%(rank))
  if True:
    #something new
    devices_per_node = int(os.environ["DEVICES_PER_NODE"])
    SIM.device_Id = rank%devices_per_node

    assert gpu_channels_singleton.get_deviceID()==SIM.device_Id
    if gpu_channels_singleton.get_nchannels() == 0: # if uninitialized
        P = Profiler("Initialize the channels singleton rank %d"%(rank))
        for x in range(len(flux)):
          gpu_channels_singleton.structure_factors_to_GPU_direct(
           x, sfall_channels[x].indices(), sfall_channels[x].data())
        del P
        import time
        print("datetime for channels singleton rank %d"%(rank),time.time())

    gpu_simulation = exascale_api(nanoBragg = SIM)
    gpu_simulation.allocate_cuda()

    gpu_detector = gpud(deviceId=SIM.device_Id, nanoBragg=SIM)
    gpu_detector.each_image_allocate()

    # loop over energies
    for x in range(len(flux)):
      P = Profiler("USE_EXASCALE_API nanoBragg Python and C++ rank %d"%(rank))

      print("USE_EXASCALE_API+++++++++++++++++++++++ Wavelength",x)
      # from channel_pixels function
      SIM.wavelength_A = wavlen[x]
      SIM.flux = flux[x]
      gpu_simulation.add_energy_channel_from_gpu_amplitudes_cuda(
        x, gpu_channels_singleton, gpu_detector)
      del P
    gpu_detector.scale_in_place(crystal.domains_per_crystal) # apply scale directly on GPU
    SIM.wavelength_A = wavelength_A # return to canonical energy for subsequent background

    if add_background_algorithm == "cuda":
      QQ = Profiler("nanoBragg background rank %d"%(rank))
      SIM.Fbg_vs_stol = water_bg
      SIM.amorphous_sample_thick_mm = 0.1
      SIM.amorphous_density_gcm3 = 1
      SIM.amorphous_molecular_weight_Da = 18
      SIM.flux=1e12
      SIM.beamsize_mm=0.003 # square (not user specified)
      SIM.exposure_s=1.0 # multiplies flux x exposure
      gpu_simulation.add_background_cuda(gpu_detector)
      SIM.Fbg_vs_stol = air_bg
      SIM.amorphous_sample_thick_mm = 10 # between beamstop and collimator
      SIM.amorphous_density_gcm3 = 1.2e-3
      SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
      gpu_simulation.add_background_cuda(gpu_detector)

    # deallocate GPU arrays
    gpu_detector.write_raw_pixels(SIM)  # updates SIM.raw_pixels from GPU
    gpu_detector.each_image_free()
    SIM.Amatrix_RUB = Amatrix_rot # return to canonical orientation
    del QQ

  if add_background_algorithm in ["jh","sort_stable"]:
    QQ = Profiler("nanoBragg background rank %d"%(rank))

    SIM.Fbg_vs_stol = water_bg
    SIM.amorphous_sample_thick_mm = 0.1
    SIM.amorphous_density_gcm3 = 1
    SIM.amorphous_molecular_weight_Da = 18
    SIM.flux=1e12
    SIM.beamsize_mm=0.003 # square (not user specified)
    SIM.exposure_s=1.0 # multiplies flux x exposure
    SIM.add_background(sort_stable=(add_background_algorithm=="sort_stable"))

    SIM.Fbg_vs_stol = air_bg
    SIM.amorphous_sample_thick_mm = 10 # between beamstop and collimator
    SIM.amorphous_density_gcm3 = 1.2e-3
    SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
    SIM.add_background(sort_stable=(add_background_algorithm=="sort_stable"))
    del QQ

  SIM.detector_psf_kernel_radius_pixels=5;
  SIM.detector_psf_type=shapetype.Unknown # for CSPAD
  SIM.detector_psf_fwhm_mm=0
  #SIM.apply_psf()

  QQ = Profiler("nanoBragg noise rank %d"%(rank))
  #SIM.add_noise() #converts phtons to ADU.
  del QQ

  extra = "PREFIX=%s;\nRANK=%d;\n"%(prefix,rank)
  SIM.to_smv_format_py(fileout=burst_buffer_fileout,intfile_scale=1,rotmat=True,extra=extra,gz=True)

  SIM.free_all()


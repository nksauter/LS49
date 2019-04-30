from __future__ import division,print_function
from six.moves import range
from six.moves import StringIO
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from cctbx import crystal_orientation
import libtbx.load_env # possibly implicit
from cctbx import crystal
import math
import scitbx
from LS49.sim.util_fmodel import gen_fmodel

big_data = "." # directory location for reference files
def full_path(filename):
  import os
  big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  return os.path.join(big_data,filename)

def data():
  from LS49.sim.fdp_plot import george_sherrell
  return dict(
    pdb_lines = open(full_path("1m2a.pdb"),"r").read(),
    Fe_oxidized_model = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out")),
    Fe_reduced_model = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out")),
    Fe_metallic_model = george_sherrell(full_path("data_sherrell/Fe_fake.dat"))
  )

from LS49.sim.step4_pad import microcrystal

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
  import os
  return (not os.path.isfile(fname)) and (not os.path.isfile(fname+".gz"))

add_spots_algorithm = "NKS"


class ChannelSimulator:
  def __init__(self, UMAT_nm, N,
               Amatrix_rot,
               mosaic_domains=25,
               mosaic_spread_deg=0.05,
               SEED=1,
               randomize=False):
    """

    :param UMAT_nm:  mosaic blocks argument
    :param N:  number of unit cells along dim
    :param Amatrix_rot: rotated Amtrix im nanoBragg format
    """
    self.SIM = nanoBragg(
                detpixels_slowfast=(3000,3000),
                pixel_size_mm=0.11,
                Ncells_abc=(N, N, N),
                wavelength_A=1,  # default we will update later@
                verbose=0)

    self.SIM.seed = SEED

    self.SIM.adc_offset_adu = 10 # Do not offset by 40
    self.SIM.mosaic_domains = mosaic_domains  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
    self.SIM.mosaic_spread_deg = mosaic_spread_deg # interpreted by UMAT_nm as a half-width stddev
    self.SIM.distance_mm=141.7
    self.SIM.distance_mm=141.7
    self.SIM.set_mosaic_blocks(UMAT_nm)
    self.SIM.oversample=1
    self.SIM.polarization=1
    self.SIM.default_F=1e5  # TODO: in the future we will init the energy dependent F_HKL here
    #self.SIM.Fhkl = energy_independent_F
    self.SIM.Amatrix_RUB = Amatrix_rot
    self.SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
    self.SIM.progress_meter=False
    self.SIM.exposure_s=1.0 # so total fluence is e12
    self.SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
    if randomize:
      self.SIM.random_orientation()
    temp=self.SIM.Ncells_abc
    print("Ncells_abc=",self.SIM.Ncells_abc)
    self.SIM.Ncells_abc=temp

    # FIXME: add the CUDA init script here
    #initialize_GPU_variables()
    self.raw_pixels = self.SIM.raw_pixels  # FIXME: this will be on GPU

  def add_channel_pixels(self, wavelength_A, flux, rank):
    """
    :param wavelength_A: wavelength float parameter in Angstrom
    :param flux: flux of photons of the given wavelength
    :param rank: job rank ?
    :return:
    """

    self.SIM.flux=flux
    self.SIM.wavelength_A = wavelength_A

    from libtbx.development.timers import Profiler
    P = Profiler("nanoBragg C++ rank %d"%(rank))
    if add_spots_algorithm is "NKS":
      from boost.python import streambuf # will deposit printout into dummy StringIO as side effect
      self.SIM.add_nanoBragg_spots_nks(streambuf(StringIO()))
    elif add_spots_algorithm is "JH":
      self.SIM.add_nanoBragg_spots()
    elif add_spots_algorithm is "cuda":
      self.SIM.add_nanoBragg_spots_cuda()  # TODO: figure out if this sets pixels to 0 or SIM.raw_pixels
      #self.SIM.add_nanoBragg_spots_cuda_light()  # FIXME: this should be updated to not initialize

    else: raise Exception("unknown spots algorithm")
    del P

    self.raw_pixels += self.SIM.raw_pixels

from LS49.sim.debug_utils import channel_extractor
CHDBG_singleton = channel_extractor()

def run_sim2smv(prefix,crystal,spectra,rotation,rank,quick=False):
  local_data = data()
  smv_fileout = prefix + ".img"
  if quick is not True:
    if not write_safe(smv_fileout):
      print("File %s already exists, skipping in rank %d"%(smv_fileout,rank))
      return

  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = next(spectra) # list of lambdas, list of fluxes, average wavelength
  assert wavelength_A > 0
  if quick:
    wavlen = flex.double([wavelength_A]);
    flux = flex.double([flex.sum(flux)])
    print("Quick sim, lambda=%f, flux=%f"%(wavelength_A,flux[0]))

  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()
  sfall_main = GF.get_amplitudes()

  # use crystal structure to initialize Fhkl array
  sfall_main.show_summary(prefix = "Amplitudes used ")
  N = crystal.number_of_cells(sfall_main.unit_cell())

  import sys

  mosaic_spread_deg = 0.05
  mosaic_domains = 25
  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  scitbx.random.set_random_seed(1234)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=mosaic_spread_deg * math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(mosaic_domains)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )

  #SIM.detector_thick_mm = 0.5 # = 0 for Rayonix
  #SIM.detector_thicksteps = 1 # should default to 1 for Rayonix, but set to 5 for CSPAD
  #SIM.detector_attenuation_length_mm = default is silicon
  Amatrix_rot = (rotation * sqr(sfall_main.unit_cell().orthogonalization_matrix())).transpose()

  # get same noise each time this test is run
  # this will become F000, marking the beam center
  #SIM.missets_deg= (10,20,30)

  channel_sim = ChannelSimulator(
    N=N,
    UMAT_nm=UMAT_nm,
    Amatrix_rot=Amatrix_rot)

  for x in range(len(flux)):

    from libtbx.development.timers import Profiler
    P = Profiler("nanoBragg Python and C++ rank %d"%(rank))

    print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
    channel_sim.add_channel_pixels(wavlen[x], flux[x], rank)
    CHDBG_singleton.extract(channel_no=x, data=channel_sim.raw_pixels)
    del P

  SIM = channel_sim.SIM

  # FIXME: this hsould already be done by kernel...
  SIM.raw_pixels = channel_sim.raw_pixels * crystal.domains_per_crystal;

  if quick:  SIM.to_smv_format(fileout=prefix + "_intimage_001.img")

  bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 0.1
  SIM.amorphous_density_gcm3 = 1
  SIM.amorphous_molecular_weight_Da = 18
  SIM.flux=1e12
  SIM.beamsize_mm=0.003 # square (not user specified)
  SIM.exposure_s=1.0 # multiplies flux x exposure
  SIM.add_background()
  if quick:  SIM.to_smv_format(fileout=prefix + "_intimage_002.img")

  # rough approximation to air
  bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),(0.5,4.22)])
  SIM.Fbg_vs_stol = bg
  #SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
  SIM.amorphous_sample_thick_mm = 10 # between beamstop and collimator
  SIM.amorphous_density_gcm3 = 1.2e-3
  SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
  print("amorphous_sample_size_mm=",SIM.amorphous_sample_size_mm)
  print("amorphous_sample_thick_mm=",SIM.amorphous_sample_thick_mm)
  print("amorphous_density_gcm3=",SIM.amorphous_density_gcm3)
  print("amorphous_molecular_weight_Da=",SIM.amorphous_molecular_weight_Da)
  SIM.add_background()

  #apply beamstop mask here

  # set this to 0 or -1 to trigger automatic radius.  could be very slow with bright images
  # settings for CCD
  SIM.detector_psf_kernel_radius_pixels=5;
  #SIM.detector_psf_fwhm_mm=0.08;
  #SIM.detector_psf_type=shapetype.Fiber # rayonix=Fiber, CSPAD=None (or small Gaussian)
  SIM.detector_psf_type=shapetype.Unknown # for CSPAD
  SIM.detector_psf_fwhm_mm=0
  #SIM.apply_psf()
  print("One pixel-->",SIM.raw_pixels[500000])

  # at this point we scale the raw pixels so that the output array is on an scale from 0 to 50000.
  # that is the default behavior (intfile_scale<=0), otherwise it applies intfile_scale as a multiplier on an abs scale.
  if quick:
    SIM.to_smv_format(fileout=prefix + "_intimage_003.img")

  print("quantum_gain=",SIM.quantum_gain) #defaults to 1. converts photons to ADU
  print("adc_offset_adu=",SIM.adc_offset_adu)
  print("detector_calibration_noise_pct=",SIM.detector_calibration_noise_pct)
  print("flicker_noise_pct=",SIM.flicker_noise_pct)
  print("readout_noise_adu=",SIM.readout_noise_adu) # gaussian random number to add to every pixel (0 for PAD)
  # apply Poissonion correction, then scale to ADU, then adc_offset.
  # should be 10 for most Rayonix, Pilatus should be 0, CSPAD should be 0.

  print("detector_psf_type=",SIM.detector_psf_type)
  print("detector_psf_fwhm_mm=",SIM.detector_psf_fwhm_mm)
  print("detector_psf_kernel_radius_pixels=",SIM.detector_psf_kernel_radius_pixels)
  SIM.add_noise() #converts phtons to ADU.

  print("raw_pixels=",SIM.raw_pixels)
  extra = "PREFIX=%s;\nRANK=%d;\n"%(prefix,rank)
  SIM.to_smv_format_py(fileout=smv_fileout,intfile_scale=1,rotmat=True,extra=extra,gz=True)

  # try to write as CBF
  if False:
    import dxtbx
    from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
    img = dxtbx.load(prefix + ".img")
    print(img)
    FormatCBFMiniPilatus.as_file(
    detector=img.get_detector(),beam=img.get_beam(),gonio=img.get_goniometer(),scan=img.get_scan(),
    data=img.get_raw_data(),path=prefix + ".cbf")
  SIM.free_all()

def tst_all(quick=False,prefix="step5"):
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_images(20,energy=7120.,total_flux=1e12)

  #
  C = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
  mt = flex.mersenne_twister(seed=0)

  if quick: prefix_root=prefix + "_%06d"
  else: prefix_root=prefix + "poly_%06d"

  Nimages = 1# 10000
  for iteration in range(Nimages):
    file_prefix = prefix_root%iteration
    rand_ori = sqr(mt.random_double_r3_rotation_matrix())
    run_sim2smv(prefix = file_prefix,crystal = C,spectra=iterator,rotation=rand_ori,quick=quick,rank=0)

if __name__=="__main__":
  tst_all()
  print("OK")

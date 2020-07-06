from __future__ import division,print_function
from six.moves import range
from scitbx.array_family import flex
from scitbx.matrix import sqr
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step4_pad import microcrystal
from libtbx.development.timers import Profiler
from LS49.sim.channel_simulator import  ChannelSimulator


# GLOBAL
ADD_SPOTS_ALGORITHM = "NKS"

def full_path(filename):
  import os
  from LS49 import ls49_big_data
  return os.path.join(ls49_big_data,filename)

def data():
  from LS49.sim.fdp_plot import george_sherrell
  return dict(
    pdb_lines = open(full_path("1m2a.pdb"),"r").read(),
    Fe_oxidized_model = george_sherrell(full_path("data_sherrell/pf-rd-ox_fftkk.out")),
    Fe_reduced_model = george_sherrell(full_path("data_sherrell/pf-rd-red_fftkk.out")),
    Fe_metallic_model = george_sherrell(full_path("data_sherrell/Fe_fake.dat"))
  )


def write_safe(fname):
  # make sure file or compressed file is not already on disk
  import os
  return (not os.path.isfile(fname)) and (not os.path.isfile(fname+".gz"))


from LS49.sim.debug_utils import channel_extractor
CHDBG_singleton = channel_extractor()

def run_sim2smv(prefix,crystal,spectra,rotation,rank):
  local_data = data()
  smv_fileout = prefix + ".img"
  if not write_safe(smv_fileout):
    print("File %s already exists, skipping in rank %d"%(smv_fileout,rank))
    return

  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = next(spectra) # list of lambdas, list of fluxes, average wavelength
  assert wavelength_A > 0

  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()
  sfall_main = GF.get_amplitudes()

  # use crystal structure to initialize Fhkl array
  sfall_main.show_summary(prefix = "Amplitudes used ")
  N = crystal.number_of_cells(sfall_main.unit_cell())
  mosaic_spread_deg = 0.05
  mosaic_domains = 25

  channel_sim = ChannelSimulator(
    rotation=rotation,
    sfall_main=sfall_main,
    N=N,
    mosaic_spread_deg=mosaic_spread_deg,
    mosaic_domains=mosaic_domains)

  for x in range(len(flux)):

    P = Profiler("nanoBragg Python and C++ rank %d"%(rank))

    print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
    channel_sim.add_channel_pixels(wavlen[x], flux[x], rank,
                                   algo=ADD_SPOTS_ALGORITHM)

    # whats this ?
    CHDBG_singleton.extract(channel_no=x, data=channel_sim.raw_pixels)
    del P

  SIM = channel_sim.SIM

  # FIXME: this hsould already be done by kernel...
  SIM.raw_pixels = channel_sim.raw_pixels*crystal.domains_per_crystal

  bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 0.1
  SIM.amorphous_density_gcm3 = 1
  SIM.amorphous_molecular_weight_Da = 18
  SIM.flux=1e12
  SIM.beamsize_mm=0.003 # square (not user specified)
  SIM.exposure_s=1.0 # multiplies flux x exposure
  SIM.add_background()

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

  SIM.free_all()

def tst_all(prefix="step5"):
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_images(20,energy=7120.,total_flux=1e12)
  #
  C = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
  mt = flex.mersenne_twister(seed=0)

  prefix_root=prefix + "poly_%06d"

  Nimages = 1# 10000
  for iteration in range(Nimages):
    file_prefix = prefix_root%iteration
    rand_ori = sqr(mt.random_double_r3_rotation_matrix())
    run_sim2smv(prefix = file_prefix,crystal = C,spectra=iterator,rotation=rand_ori,rank=0)

if __name__=="__main__":
  tst_all()
  print("OK")

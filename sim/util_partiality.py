from __future__ import division, print_function
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.array_family import flex
from LS49.sim.step5_pad import microcrystal
from scitbx.matrix import sqr,col
from LS49.sim.step5_pad import pdb_lines
from LS49.sim.util_fmodel import gen_fmodel
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from six.moves import StringIO
import scitbx
import math
from LS49.sim.fdp_plot import george_sherrell
Fe_oxidized_model = george_sherrell("data_sherrell/pf-rd-ox_fftkk.out")
Fe_reduced_model = george_sherrell("data_sherrell/pf-rd-red_fftkk.out")
global g_sfall
g_sfall = None
use_g_sfall=False

# use local file with (open("sfall_P1_7122_amplitudes.pickle","wb")) as F:
with (open("sfall_P1_7122_amplitudes.pickle","rb")) as F:
  sfall_7122 = pickle.load(F)

def channel_pixels(ROI,wavelength_A,flux,N,UMAT_nm,Amatrix_rot,fmodel_generator,output):
  energy_dependent_fmodel=False
  if energy_dependent_fmodel:
    fmodel_generator.reset_wavelength(wavelength_A)
    fmodel_generator.reset_specific_at_wavelength(
                   label_has="FE1",tables=Fe_oxidized_model,newvalue=wavelength_A)
    fmodel_generator.reset_specific_at_wavelength(
                   label_has="FE2",tables=Fe_reduced_model,newvalue=wavelength_A)
    print("USING scatterer-specific energy-dependent scattering factors")
    sfall_channel = fmodel_generator.get_amplitudes()
  elif use_g_sfall:
    global g_sfall
    if g_sfall is None:  g_sfall = fmodel_generator.get_amplitudes()
    sfall_channel = g_sfall
  else:
    sfall_channel = sfall_7122

  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 10 # Do not offset by 40
  SIM.region_of_interest = ROI
  SIM.mosaic_domains = 50  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
  SIM.distance_mm=141.7
  SIM.set_mosaic_blocks(UMAT_nm)

  #SIM.detector_thick_mm = 0.5 # = 0 for Rayonix
  #SIM.detector_thicksteps = 1 # should default to 1 for Rayonix, but set to 5 for CSPAD
  #SIM.detector_attenuation_length_mm = default is silicon

  # get same noise each time this test is run
  SIM.seed = 1
  SIM.oversample=1
  SIM.wavelength_A = wavelength_A
  SIM.polarization=1
  SIM.default_F=0
  SIM.Fhkl=sfall_channel
  SIM.Amatrix_RUB = Amatrix_rot
  SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
  SIM.progress_meter=False
  # flux is always in photons/s
  SIM.flux=flux
  SIM.exposure_s=1.0 # so total fluence is e12
  # assumes round beam
  SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
  temp=SIM.Ncells_abc
  print("Ncells_abc=",SIM.Ncells_abc)
  SIM.Ncells_abc=temp

  from libtbx.development.timers import Profiler
  P = Profiler("nanoBragg")
  SIM.printout=True
  fast = ROI[1][0] + (ROI[1][1]-ROI[1][0])//2
  slow = ROI[0][0] + (ROI[0][1]-ROI[0][0])//2
  SIM.printout_pixel_fastslow=(slow,fast)
  from boost.python import streambuf # will deposit printout into output as side effect
  SIM.add_nanoBragg_spots_nks(streambuf(output))
  del P
  print(("SIM count > 0",(SIM.raw_pixels>0).count(True)))
  return SIM

def run_sim2smv(ROI,prefix,crystal,spectra,rotation,rank,quick=False):
  smv_fileout = prefix + ".img"

  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = next(spectra) # list of lambdas, list of fluxes, average wavelength

  tophat_spectrum = True
  if tophat_spectrum:
    sum_flux = flex.sum(flux)
    #from IPython import embed; embed()
    ave_flux = sum_flux/60. # 60 energy channels
    for ix in range(len(wavlen)):
      energy = 12398.425 / wavlen[ix]
      if energy>=7090 and energy <=7150:
        flux[ix]=ave_flux
      else:
        flux[ix]=0.
  if quick:
    wavlen = flex.double([wavelength_A]);
    flux = flex.double([flex.sum(flux)])
    print("Quick sim, lambda=%f, flux=%f"%(wavelength_A,flux[0]))

  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=pdb_lines,algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()
  sfall_main = GF.get_amplitudes()

  # use crystal structure to initialize Fhkl array
  sfall_main.show_summary(prefix = "Amplitudes used ")
  N = crystal.number_of_cells(sfall_main.unit_cell())

  #SIM = nanoBragg(detpixels_slowfast=(2000,2000),pixel_size_mm=0.11,Ncells_abc=(5,5,5),verbose=0)
  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    # workaround for problem with wavelength array, specify it separately in constructor.
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 0 # Do not offset by 40
  SIM.adc_offset_adu = 10 # Do not offset by 40
  import sys
  if len(sys.argv)>2:
    SIM.seed = -int(sys.argv[2])
    print("GOTHERE seed=",SIM.seed)
  if len(sys.argv)>1:
    if sys.argv[1]=="random" : SIM.randomize_orientation()
  SIM.mosaic_domains = 50  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
                           # 3000000 images would be 100000 hours on a 60-core machine (dials), or 11.4 years
                           # using 2 nodes, 5.7 years.  Do this at SLAC? NERSC? combination of all?
                           # SLAC downtimes: Tues Dec 5 (24 hrs), Mon Dec 11 (72 hrs), Mon Dec 18 light use, 24 days
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
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
  #SIM.missets_deg= (10,20,30)
  print("mosaic_seed=",SIM.mosaic_seed)
  print("seed=",SIM.seed)
  print("calib_seed=",SIM.calib_seed)
  print("missets_deg =", SIM.missets_deg)
  SIM.Fhkl=sfall_main
  print("Determinant",rotation.determinant())
  Amatrix_rot = (rotation * sqr(sfall_main.unit_cell().orthogonalization_matrix())).transpose()
  print("RAND_ORI", prefix, end=' ')
  for i in Amatrix_rot: print(i, end=' ')
  print()

  SIM.Amatrix_RUB = Amatrix_rot
  #workaround for failing init_cell, use custom written Amatrix setter
  print("unit_cell_Adeg=",SIM.unit_cell_Adeg)
  print("unit_cell_tuple=",SIM.unit_cell_tuple)
  Amat = sqr(SIM.Amatrix).transpose() # recovered Amatrix from SIM
  from cctbx import crystal_orientation
  Ori = crystal_orientation.crystal_orientation(Amat, crystal_orientation.basis_type.reciprocal)
  print("Python unit cell from SIM state",Ori.unit_cell())

  # fastest option, least realistic
  #SIM.xtal_shape=shapetype.Tophat # RLP = hard sphere
  #SIM.xtal_shape=shapetype.Square # gives fringes
  SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
  #SIM.xtal_shape=shapetype.Round # Crystal is a hard sphere
  # only really useful for long runs
  SIM.progress_meter=False
  # prints out value of one pixel only.  will not render full image!
  #SIM.printout_pixel_fastslow=(500,500)
  #SIM.printout=True
  SIM.show_params()
  # flux is always in photons/s
  SIM.flux=1e12
  SIM.exposure_s=1.0 # so total fluence is e12
  # assumes round beam
  SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
  temp=SIM.Ncells_abc
  print("Ncells_abc=",SIM.Ncells_abc)
  SIM.Ncells_abc=temp
  print("Ncells_abc=",SIM.Ncells_abc)
  print("xtal_size_mm=",SIM.xtal_size_mm)
  print("unit_cell_Adeg=",SIM.unit_cell_Adeg)
  print("unit_cell_tuple=",SIM.unit_cell_tuple)
  print("missets_deg=",SIM.missets_deg)
  print("Amatrix=",SIM.Amatrix)
  print("beam_center_mm=",SIM.beam_center_mm)
  print("XDS_ORGXY=",SIM.XDS_ORGXY)
  print("detector_pivot=",SIM.detector_pivot)
  print("xtal_shape=",SIM.xtal_shape)
  print("beamcenter_convention=",SIM.beamcenter_convention)
  print("fdet_vector=",SIM.fdet_vector)
  print("sdet_vector=",SIM.sdet_vector)
  print("odet_vector=",SIM.odet_vector)
  print("beam_vector=",SIM.beam_vector)
  print("polar_vector=",SIM.polar_vector)
  print("spindle_axis=",SIM.spindle_axis)
  print("twotheta_axis=",SIM.twotheta_axis)
  print("distance_meters=",SIM.distance_meters)
  print("distance_mm=",SIM.distance_mm)
  print("close_distance_mm=",SIM.close_distance_mm)
  print("detector_twotheta_deg=",SIM.detector_twotheta_deg)
  print("detsize_fastslow_mm=",SIM.detsize_fastslow_mm)
  print("detpixels_fastslow=",SIM.detpixels_fastslow)
  print("detector_rot_deg=",SIM.detector_rot_deg)
  print("curved_detector=",SIM.curved_detector)
  print("pixel_size_mm=",SIM.pixel_size_mm)
  print("point_pixel=",SIM.point_pixel)
  print("polarization=",SIM.polarization)
  print("nopolar=",SIM.nopolar)
  print("oversample=",SIM.oversample)
  print("region_of_interest=",SIM.region_of_interest)
  print("wavelength_A=",SIM.wavelength_A)
  print("energy_eV=",SIM.energy_eV)
  print("fluence=",SIM.fluence)
  print("flux=",SIM.flux)
  print("exposure_s=",SIM.exposure_s)
  print("beamsize_mm=",SIM.beamsize_mm)
  print("dispersion_pct=",SIM.dispersion_pct)
  print("dispsteps=",SIM.dispsteps)
  print("divergence_hv_mrad=",SIM.divergence_hv_mrad)
  print("divsteps_hv=",SIM.divsteps_hv)
  print("divstep_hv_mrad=",SIM.divstep_hv_mrad)
  print("round_div=",SIM.round_div)
  print("phi_deg=",SIM.phi_deg)
  print("osc_deg=",SIM.osc_deg)
  print("phisteps=",SIM.phisteps)
  print("phistep_deg=",SIM.phistep_deg)
  print("detector_thick_mm=",SIM.detector_thick_mm)
  print("detector_thicksteps=",SIM.detector_thicksteps)
  print("detector_thickstep_mm=",SIM.detector_thickstep_mm)
  print("***mosaic_spread_deg=",SIM.mosaic_spread_deg)
  print("***mosaic_domains=",SIM.mosaic_domains)
  print("indices=",SIM.indices)
  print("amplitudes=",SIM.amplitudes)
  print("Fhkl_tuple=",SIM.Fhkl_tuple)
  print("default_F=",SIM.default_F)
  print("interpolate=",SIM.interpolate)
  print("integral_form=",SIM.integral_form)

  from libtbx.development.timers import Profiler
  P = Profiler("nanoBragg")
  # now actually burn up some CPU
  #SIM.add_nanoBragg_spots()
  del P

  # simulated crystal is only 125 unit cells (25 nm wide)
  # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
  print(crystal.domains_per_crystal)
  SIM.raw_pixels *= crystal.domains_per_crystal; # must calculate the correct scale!
  output = StringIO() # open("myfile","w")
  for x in range(0,100,2): #len(flux)):
    if flux[x]==0.0:continue
    print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
    CH = channel_pixels(ROI,wavlen[x],flux[x],N,UMAT_nm,Amatrix_rot,GF,output)
    SIM.raw_pixels += CH.raw_pixels * crystal.domains_per_crystal;
    print(SIM.raw_pixels)

    CH.free_all()

  message = output.getvalue().split()
  miller = (int(message[4]),int(message[5]),int(message[6]))
  intensity = float(message[9]);

  #SIM.to_smv_format(fileout=prefix + "_intimage_001.img")
  pixels = SIM.raw_pixels
  roi_pixels = pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
  print("Reducing full shape of",pixels.focus(),"to ROI of",roi_pixels.focus())
  SIM.free_all()
  return dict(roi_pixels=roi_pixels,miller=miller,intensity=intensity)

def get_partiality_response(key,one_index,spectra_simulation,ROI):

  N_total = 100000 # number of items to simulate
  spectra = spectra_simulation
  crystal = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
  mt = flex.mersenne_twister(seed=0)
  random_orientations = []
  for iteration in range(N_total):
    random_orientations.append( mt.random_double_r3_rotation_matrix() )

  iterator = spectra.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)

  file_prefix = "key_slow_nonoise_%06d"%key
  rand_ori = sqr(random_orientations[key])

  pixels = run_sim2smv(ROI,prefix = file_prefix,crystal = crystal,spectra=iterator,rotation=rand_ori,quick=False,rank=0)
  return pixels

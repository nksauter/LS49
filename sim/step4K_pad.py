from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal
import math
import scitbx
from LS49.sim.util_fmodel import fmodel_from_pdb, fcalc_from_pdb

pdb_lines = open("1m2a.pdb","r").read()
#pdb_lines = open("4tnl.pdb","r").read()

from LS49.sim.step4_pad import microcrystal

"""Changes in Step4K relative to Step4
k_sol is now hard-coded to 0.435 instead of 0.35 to reduce solvent contrast
Simulation out to 1.7 Angstrom resolution
Use 3K x 3K detector, not 2K x 2K to insure complete coverage at 2.1 A inscribed circle
Put the ROTMAT orientation in the file header, instead of in log file
Put sequence number and MPI rank in header
Write gz compressed file
Avoid computing the image if already on disk
"""

"""Prospective plan
Simulate an anomalous dataset, 7150 eV, but using fixed f',f". 50,000 images. Solve by SAD. (step4K)
*** Plug in the f' & f" as a function of wavelength.  Rerun the simulation at 7120 eV.
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

def channel_pixels(wavelength_A,flux,N,UMAT_nm,Amatrix_rot,sfall):
  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 10 # Do not offset by 40
  SIM.mosaic_domains = 25  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
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
  SIM.Fhkl=sfall
  SIM.Amatrix_RUB = Amatrix_rot
  SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
  SIM.progress_meter=False
  # flux is always in photons/s
  SIM.flux=flux
  SIM.exposure_s=1.0 # so total fluence is e12
  # assumes round beam
  SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
  temp=SIM.Ncells_abc
  print "Ncells_abc=",SIM.Ncells_abc
  SIM.Ncells_abc=temp

  from libtbx.development.timers import Profiler
  P = Profiler("nanoBragg")
  SIM.add_nanoBragg_spots_nks()
  del P
  return SIM

def run_sim2smv(prefix,crystal,spectra,rotation,rank,quick=False):
  smv_fileout = prefix + ".img"
  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = spectra.next() # list of lambdas, list of fluxes, average wavelength
  if quick:
    wavlen = flex.double([wavelength_A]);
    flux = flex.double([flex.sum(flux)])
    print "Quick sim, lambda=%f, flux=%f"%(wavelength_A,flux[0])

  #sfall = fcalc_from_pdb(resolution=direct_algo_res_limit,pdb_text=pdb_lines,algorithm="direct",wavelength=SIM.wavelength_A)
  sfall = fmodel_from_pdb(resolution=direct_algo_res_limit,pdb_text=pdb_lines,algorithm="fft",wavelength=wavelength_A)

  # use crystal structure to initialize Fhkl array
  sfall.show_summary(prefix = "Amplitudes used ")
  N = crystal.number_of_cells(sfall.unit_cell())

  #SIM = nanoBragg(detpixels_slowfast=(2000,2000),pixel_size_mm=0.11,Ncells_abc=(5,5,5),verbose=0)
  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    # workaround for problem with wavelength array, specify it separately in constructor.
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 0 # Do not offset by 40
  SIM.adc_offset_adu = 10 # Do not offset by 40
  import sys
  if len(sys.argv)>2:
    SIM.seed = -int(sys.argv[2])
    print "GOTHERE seed=",SIM.seed
  if len(sys.argv)>1:
    if sys.argv[1]=="random" : SIM.randomize_orientation()
  SIM.mosaic_domains = 25  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
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
  print "mosaic_seed=",SIM.mosaic_seed
  print "seed=",SIM.seed
  print "calib_seed=",SIM.calib_seed
  print "missets_deg =", SIM.missets_deg
  SIM.Fhkl=sfall
  print "Determinant",rotation.determinant()
  Amatrix_rot = (rotation * sqr(sfall.unit_cell().orthogonalization_matrix())).transpose()
  print "RAND_ORI", prefix,
  for i in Amatrix_rot: print i,
  print

  SIM.Amatrix_RUB = Amatrix_rot
  #workaround for failing init_cell, use custom written Amatrix setter
  print "unit_cell_Adeg=",SIM.unit_cell_Adeg
  print "unit_cell_tuple=",SIM.unit_cell_tuple
  Amat = sqr(SIM.Amatrix).transpose() # recovered Amatrix from SIM
  from cctbx import crystal_orientation
  Ori = crystal_orientation.crystal_orientation(Amat, crystal_orientation.basis_type.reciprocal)
  print "Python unit cell from SIM state",Ori.unit_cell()

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
  print "Ncells_abc=",SIM.Ncells_abc
  SIM.Ncells_abc=temp
  print "Ncells_abc=",SIM.Ncells_abc
  print "xtal_size_mm=",SIM.xtal_size_mm
  print "unit_cell_Adeg=",SIM.unit_cell_Adeg
  print "unit_cell_tuple=",SIM.unit_cell_tuple
  print "missets_deg=",SIM.missets_deg
  print "Amatrix=",SIM.Amatrix
  print "beam_center_mm=",SIM.beam_center_mm
  print "XDS_ORGXY=",SIM.XDS_ORGXY
  print "detector_pivot=",SIM.detector_pivot
  print "xtal_shape=",SIM.xtal_shape
  print "beamcenter_convention=",SIM.beamcenter_convention
  print "fdet_vector=",SIM.fdet_vector
  print "sdet_vector=",SIM.sdet_vector
  print "odet_vector=",SIM.odet_vector
  print "beam_vector=",SIM.beam_vector
  print "polar_vector=",SIM.polar_vector
  print "spindle_axis=",SIM.spindle_axis
  print "twotheta_axis=",SIM.twotheta_axis
  print "distance_meters=",SIM.distance_meters
  print "distance_mm=",SIM.distance_mm
  print "close_distance_mm=",SIM.close_distance_mm
  print "detector_twotheta_deg=",SIM.detector_twotheta_deg
  print "detsize_fastslow_mm=",SIM.detsize_fastslow_mm
  print "detpixels_fastslow=",SIM.detpixels_fastslow
  print "detector_rot_deg=",SIM.detector_rot_deg
  print "curved_detector=",SIM.curved_detector
  print "pixel_size_mm=",SIM.pixel_size_mm
  print "point_pixel=",SIM.point_pixel
  print "polarization=",SIM.polarization
  print "nopolar=",SIM.nopolar
  print "oversample=",SIM.oversample
  print "region_of_interest=",SIM.region_of_interest
  print "wavelength_A=",SIM.wavelength_A
  print "energy_eV=",SIM.energy_eV
  print "fluence=",SIM.fluence
  print "flux=",SIM.flux
  print "exposure_s=",SIM.exposure_s
  print "beamsize_mm=",SIM.beamsize_mm
  print "dispersion_pct=",SIM.dispersion_pct
  print "dispsteps=",SIM.dispsteps
  print "divergence_hv_mrad=",SIM.divergence_hv_mrad
  print "divsteps_hv=",SIM.divsteps_hv
  print "divstep_hv_mrad=",SIM.divstep_hv_mrad
  print "round_div=",SIM.round_div
  print "phi_deg=",SIM.phi_deg
  print "osc_deg=",SIM.osc_deg
  print "phisteps=",SIM.phisteps
  print "phistep_deg=",SIM.phistep_deg
  print "detector_thick_mm=",SIM.detector_thick_mm
  print "detector_thicksteps=",SIM.detector_thicksteps
  print "detector_thickstep_mm=",SIM.detector_thickstep_mm
  print "***mosaic_spread_deg=",SIM.mosaic_spread_deg
  print "***mosaic_domains=",SIM.mosaic_domains
  print "indices=",SIM.indices
  print "amplitudes=",SIM.amplitudes
  print "Fhkl_tuple=",SIM.Fhkl_tuple
  print "default_F=",SIM.default_F
  print "interpolate=",SIM.interpolate
  print "integral_form=",SIM.integral_form

  from libtbx.development.timers import Profiler
  P = Profiler("nanoBragg")
  # now actually burn up some CPU
  #SIM.add_nanoBragg_spots()
  del P

  # simulated crystal is only 125 unit cells (25 nm wide)
  # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
  print crystal.domains_per_crystal
  SIM.raw_pixels *= crystal.domains_per_crystal; # must calculate the correct scale!

  for x in xrange(len(flux)):
    print "+++++++++++++++++++++++++++++++++++++++ Wavelength",x
    CH = channel_pixels(wavlen[x],flux[x],N,UMAT_nm,Amatrix_rot,sfall)
    SIM.raw_pixels += CH.raw_pixels * crystal.domains_per_crystal;
    CH.free_all()
  if quick:  SIM.to_smv_format(fileout=prefix + "_intimage_001.img")

  # rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
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
  print "amorphous_sample_size_mm=",SIM.amorphous_sample_size_mm
  print "amorphous_sample_thick_mm=",SIM.amorphous_sample_thick_mm
  print "amorphous_density_gcm3=",SIM.amorphous_density_gcm3
  print "amorphous_molecular_weight_Da=",SIM.amorphous_molecular_weight_Da
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
  print "One pixel-->",SIM.raw_pixels[500000]

  # at this point we scale the raw pixels so that the output array is on an scale from 0 to 50000.
  # that is the default behavior (intfile_scale<=0), otherwise it applies intfile_scale as a multiplier on an abs scale.
  if quick:  SIM.to_smv_format(fileout=prefix + "_intimage_003.img")

  print "quantum_gain=",SIM.quantum_gain #defaults to 1. converts photons to ADU
  print "adc_offset_adu=",SIM.adc_offset_adu
  print "detector_calibration_noise_pct=",SIM.detector_calibration_noise_pct
  print "flicker_noise_pct=",SIM.flicker_noise_pct
  print "readout_noise_adu=",SIM.readout_noise_adu # gaussian random number to add to every pixel (0 for PAD)
  # apply Poissonion correction, then scale to ADU, then adc_offset.
  # should be 10 for most Rayonix, Pilatus should be 0, CSPAD should be 0.

  print "detector_psf_type=",SIM.detector_psf_type
  print "detector_psf_fwhm_mm=",SIM.detector_psf_fwhm_mm
  print "detector_psf_kernel_radius_pixels=",SIM.detector_psf_kernel_radius_pixels
  SIM.add_noise() #converts phtons to ADU.

  print "raw_pixels=",SIM.raw_pixels
  extra = "PREFIX=%s;\nRANK=%d;\n"%(prefix,rank)
  SIM.to_smv_format_py(fileout=smv_fileout,intfile_scale=1,rotmat=True,extra=extra,gz=True)

  # try to write as CBF
  if False:
    import dxtbx
    from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
    img = dxtbx.load(prefix + ".img")
    print img
    FormatCBFMiniPilatus.as_file(
    detector=img.get_detector(),beam=img.get_beam(),gonio=img.get_goniometer(),scan=img.get_scan(),
    data=img.get_raw_data(),path=prefix + ".cbf")
  SIM.free_all()


def tst_all():
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  #SS.plot_recast_images(20,energy=7150.)
  iterator = SS.generate_recast_renormalized_images(20,energy=7150.,total_flux=1e12)

  #
  C = microcrystal(Deff_A = 4000, length_um = 1., beam_diameter_um = 1.0) # assume smaller than 10 um crystals
  mt = flex.mersenne_twister(seed=0)

  quick = True
  if quick: prefix_root="step4K_%06d"
  else: prefix_root="step4Kpoly_%06d"

  Nimages = 1# 10000
  for iteration in xrange(Nimages):
    file_prefix = prefix_root%iteration
    rand_ori = sqr(mt.random_double_r3_rotation_matrix())
    run_sim2smv(prefix = file_prefix,crystal = C,spectra=iterator,rotation=rand_ori,quick=quick,rank=0)

if __name__=="__main__":
  tst_all()
  print "OK"

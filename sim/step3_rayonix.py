from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import testuple
from simtbx.nanoBragg import pivot
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import convention
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
from cctbx import crystal
from cctbx import miller
import math

pdb_lines = """HEADER TEST
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""

"""Dev path to get step 3
monochromatic at first, but use mosaicity.
get the Ncell_abc size to be derived from the unit cell parameters and domain size spherical diameter (3000 Ang)
put in Fcalc from PDB
assemble polychromatic spectrum
deal with absolute levels
"""
"""
problems: pixel values -40
spot intensities flat:  (consequence of Tophat)
what is the crystal shape? Tophat vs Gauss vs Square makes a huge difference
water seemingly only diminishes bragg spots out to 2.76 Angstroms (consequence of direct transform res limit)
but it gives grey field throughout
how to specify a 500 um PAD
how to specity Rayonix
apply_noise() seems to have added the psf even when add_psf() is off! --confirmed: observed when detector_psf_fwhm_mm!=0
There seem to be hot pixels in the noisified image!
Why is the 'water ring' at 5.7 Angstroms?
"""

"""Initial try to produce simulations relevant to LS49.
step3:
produce a single PAD image
Ferredoxin, use Fcalc from PDB entry 1M2A, wild type ferredoxin from Aquifex aeolicus, 1.5 Angstrom
No anomalous signal at first.
Model a square PAD detector centered on direct beam, 2K x 2K, 0.11 mm pixels
Distance 141.7 mm, 500 um thick silicon.
Mosaicity angle is 0.05 stddev  
500 mosaic domains
Domain size is 3000Angstrom (but it doesn't figure in to calculation).
flux is average 10^12 photons per shot, but is scaled by per-shot intensity
Illuminated volume is 3x3circlex10um long.
Choose a standard orientation at first, then a random orientation.
For the first shot, choose 7150 as our target beam energy.
Break it up into 100 spectral channels
Altogether this is 500 x 100 = 50000 mono-perfect subimages.
For the spectral dispersion, use shot #1, run 209 LG36.
Use a helium atmosphere, unattenuated beam passes through.
Work out mosaic rotation ensemble.
"""

"""Prospective plan
Simulate a single std-setting PAD image, using absolute units.
Simulate a whole dataset, 20000 images.  Solve by molecular replacement
Simulate an anomalous dataset, 7150 eV, but using fixed f',f". 100,000 images.  Solve by SAD.
*** Plug in the f' & f" as a function of wavelength.  Rerun the simulation at 7120 eV.
  Sort images into normalized-intensity wavelength bins and calculate maps
*** Put a detector at far distance, and one nearby, process both simultaneously.
  How to modify nanoBragg so it takes an offset Jungfrau, CSPAD, or Rayonix
  Make sure we can index both at the same time.
  Observe the energy positioning of Bragg spots, as function of deltapsi.
*** simulate polarized spectroscopy
"""

def fcalc_from_pdb(resolution,algorithm=None,wavelength=0.9):
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None,lines = pdb_lines)
  xray_structure = pdb_inp.xray_structure_simple()
  #
  # take a detour to insist on calculating anomalous contribution of every atom
  scatterers = xray_structure.scatterers()
  for sc in scatterers:
    from cctbx.eltbx import sasaki, henke
    #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()
  # how do we do bulk solvent?
  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  fcalc = P1_primitive_xray_structure.structure_factors(
    d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
  print fcalc
  return fcalc.amplitudes()

def run_sim2smv(fileout):
  #SIM = nanoBragg(detpixels_slowfast=(2000,2000),pixel_size_mm=0.11,Ncells_abc=(5,5,5),verbose=0)
  SIM = nanoBragg(detpixels_slowfast=(2000,2000),pixel_size_mm=0.11,Ncells_abc=(20,20,20),verbose=0)
  SIM.adc_offset_adu = 0 # Do not offset by 40
  import sys
  if len(sys.argv)>2:
    SIM.seed = -int(sys.argv[2])
    print "GOTHERE seed=",SIM.seed
  if len(sys.argv)>1:
    if sys.argv[1]=="random" : SIM.randomize_orientation()
  SIM.mosaic_domains = 100 # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
                           # 3000000 images would be 100000 hours on a 60-core machine (dials), or 11.4 years
                           # using 2 nodes, 5.7 years.  Do this at SLAC? NERSC? combination of all?
                           # SLAC downtimes: Tues Dec 5 (24 hrs), Mon Dec 11 (72 hrs), Mon Dec 18 light use, 24 days
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
  SIM.distance_mm=141.7

  import scitbx
  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=SIM.mosaic_spread_deg * math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(SIM.mosaic_domains)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )
  SIM.set_mosaic_blocks(UMAT_nm)

  # get same noise each time this test is run
  SIM.seed = 1
  SIM.oversample=1
  SIM.wavelength_A=1.734045 # 7150 eV
  SIM.polarization=1
  #SIM.unit_cell_tuple=(50,50,50,90,90,90)
  print "unit_cell_Adeg=",SIM.unit_cell_Adeg
  print "unit_cell_tuple=",SIM.unit_cell_tuple
  # this will become F000, marking the beam center
  SIM.default_F=100
  #SIM.missets_deg= (10,20,30)
  print "mosaic_seed=",SIM.mosaic_seed
  print "seed=",SIM.seed
  print "calib_seed=",SIM.calib_seed
  print "missets_deg =", SIM.missets_deg
  direct_algo_res_limit = 1.0 #critical here to extend far enough
  sfall = fcalc_from_pdb(resolution=direct_algo_res_limit,algorithm="direct",wavelength=SIM.wavelength_A)
  # use crystal structure to initialize Fhkl array
  SIM.Fhkl=sfall
  # fastest option, least realistic
  #SIM.xtal_shape=shapetype.Tophat
  #SIM.xtal_shape=shapetype.Square # gives fringes
  SIM.xtal_shape=shapetype.Gauss
  # only really useful for long runs
  SIM.progress_meter=False
  # prints out value of one pixel only.  will not render full image!
  #SIM.printout_pixel_fastslow=(500,500)
  #SIM.printout=True
  SIM.show_params()
  # flux is always in photons/s
  SIM.flux=1e12
  # assumes round beam
  SIM.beamsize_mm=0.1
  SIM.exposure_s=0.1
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
  SIM.add_nanoBragg_spots()
  del P

  # simulated crystal is only 125 unit cells (25 nm wide)
  # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
  SIM.raw_pixels *= 64e9;
  SIM.to_smv_format(fileout="step3image_001.img")
  
  # rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
  bg = flex.vec2_double([(0,2.57),(0.0365,2.58),(0.07,2.8),(0.12,5),(0.162,8),(0.2,6.75),(0.18,7.32),(0.216,6.75),(0.236,6.5),(0.28,4.5),(0.3,4.3),(0.345,4.36),(0.436,3.77),(0.5,3.17)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 0.1
  SIM.amorphous_density_gcm3 = 1
  SIM.amorphous_molecular_weight_Da = 18
  SIM.flux=1e12
  SIM.beamsize_mm=0.1
  SIM.exposure_s=0.1
  SIM.add_background()
  SIM.to_smv_format(fileout="step3image_002.img")

  # rough approximation to air
  bg = flex.vec2_double([(0,14.1),(0.045,13.5),(0.174,8.35),(0.35,4.78),(0.5,4.22)])
  SIM.Fbg_vs_stol = bg
  SIM.amorphous_sample_thick_mm = 35 # between beamstop and collimator
  SIM.amorphous_density_gcm3 = 1.2e-3
  SIM.amorphous_sample_molecular_weight_Da = 28 # nitrogen = N2
  print "amorphous_sample_size_mm=",SIM.amorphous_sample_size_mm
  print "amorphous_sample_thick_mm=",SIM.amorphous_sample_thick_mm
  print "amorphous_density_gcm3=",SIM.amorphous_density_gcm3
  print "amorphous_molecular_weight_Da=",SIM.amorphous_molecular_weight_Da
  SIM.add_background()
  # set this to 0 or -1 to trigger automatic radius.  could be very slow with bright images
  # settings for CCD
  SIM.detector_psf_kernel_radius_pixels=5;
  #SIM.detector_psf_fwhm_mm=0.08;
  SIM.detector_psf_type=shapetype.Fiber
  SIM.detector_psf_fwhm_mm=0
  #SIM.apply_psf()
  print "One pixel-->",SIM.raw_pixels[500000]
  SIM.to_smv_format(fileout="step3image_003.img")

  print "quantum_gain=",SIM.quantum_gain
  print "adc_offset_adu=",SIM.adc_offset_adu
  print "detector_calibration_noise_pct=",SIM.detector_calibration_noise_pct
  print "flicker_noise_pct=",SIM.flicker_noise_pct
  print "readout_noise_adu=",SIM.readout_noise_adu
  print "detector_psf_type=",SIM.detector_psf_type
  print "detector_psf_fwhm_mm=",SIM.detector_psf_fwhm_mm
  print "detector_psf_kernel_radius_pixels=",SIM.detector_psf_kernel_radius_pixels
  SIM.add_noise()
  #fileout = "intimage_001.img"
  print "raw_pixels=",SIM.raw_pixels
  SIM.to_smv_format(fileout="step3noiseimage_001.img",intfile_scale=1)

  # try to write as CBF
  import dxtbx
  from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
  img = dxtbx.load("step3noiseimage_001.img")
  print img
  FormatCBFMiniPilatus.as_file(
    detector=img.get_detector(),beam=img.get_beam(),gonio=img.get_goniometer(),scan=img.get_scan(),
    data=img.get_raw_data(),path=fileout)
  SIM.free_all()


def tst_all():
  F = testuple()
  assert F == (1,2,3,4)
  #
  fileout = "step3noiseimage_001.cbf"
  #
  run_sim2smv(fileout)
  import os
  assert os.path.isfile(fileout)

  #simulation is complete, now we'll autoindex the image fragment and verify
  # that the indexed cell is similar to the input cell.


if __name__=="__main__":
  tst_all()
  print "OK"

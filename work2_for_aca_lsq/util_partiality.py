from __future__ import division, print_function
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.array_family import flex
from LS49.sim.step5_pad import microcrystal
from scitbx.matrix import sqr,col
from LS49.sim.step5_pad import data
from LS49.sim.util_fmodel import gen_fmodel
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from six.moves import StringIO
from cctbx import crystal_orientation
import scitbx
import math

import os
model_mode = os.environ["MODEL_MODE"]
# "superpower_postrefine" | "dials_refine" | "coarse_ground_truth"


# use local file with (open(something,"wb")) as F:
with (open("confirm_sfall_P1_7122_amplitudes.pickle","rb")) as F:
  sfall_P1_7122_amplitudes = pickle.load(F)

def channel_pixels(ROI,wavelength_A,flux,N,UMAT_nm,Amatrix_rot,fmodel_generator,output):
  local_data = data()
  energy_dependent_fmodel=False
  if energy_dependent_fmodel:
    fmodel_generator.reset_wavelength(wavelength_A)
    fmodel_generator.reset_specific_at_wavelength(
                   label_has="FE1",tables=local_data.get("Fe_oxidized_model"),newvalue=wavelength_A)
    fmodel_generator.reset_specific_at_wavelength(
                   label_has="FE2",tables=local_data.get("Fe_reduced_model"),newvalue=wavelength_A)
    print("USING scatterer-specific energy-dependent scattering factors")
    sfall_channel = fmodel_generator.get_amplitudes()
  else:
    sfall_channel = sfall_P1_7122_amplitudes

  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
    wavelength_A=wavelength_A,verbose=0)
  SIM.adc_offset_adu = 10 # Do not offset by 40
  SIM.region_of_interest = ROI
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
  SIM.mosaic_domains = 50  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
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

def ersatz_postrefined(idx,CB_OP_C_P,old_Amat):
  Amat = old_Amat
  from LS49.work_pre_experiment.post5_ang_misset import parse_postrefine # from teh postrefinement log "ASTAR"
  R = parse_postrefine()
  print ("get ersatz_postrefined with index",idx)
  C = crystal_orientation.crystal_orientation(Amat,crystal_orientation.basis_type.direct)
  assert R.get(idx,None) is not None
    # work up the crystal model from postrefinement
  direct_A = R[idx].inverse()
  permute = sqr((0,0,1,0,1,0,-1,0,0))
  sim_compatible = direct_A*permute # permute columns when post multiplying
  P = crystal_orientation.crystal_orientation(sim_compatible,crystal_orientation.basis_type.direct)
  C.show(legend="ground truth, P1")
  P.show(legend="postrefined, C2")
  PR = P.change_basis(CB_OP_C_P)
  PR.show(legend="Postrefined, primitive setting")

  cb_op_align = PR.best_similarity_transformation(C,200,1)
  align_PR = PR.change_basis(sqr(cb_op_align))
  align_PR.show(legend="postrefined, aligned")
  print("alignment matrix", cb_op_align)
  return align_PR.direct_matrix()

def superpower_postrefine(idx,CB_OP_C_P,old_Amat):
  Amat = old_Amat
  from LS49.work_pre_experiment.post5_ang_misset import get_item
  from LS49.ML_push.exploratory_missetting import metric
  R = get_item(idx)
  print ("get ersatz_postrefined with index",idx)
  C = crystal_orientation.crystal_orientation(
      Amat,crystal_orientation.basis_type.direct)
  C.show(legend="ground truth, P1")
  C2 = C.change_basis(CB_OP_C_P.inverse())
  C2.show(legend="ground truth, C2")
  direct_A = R["integrated_crystal_model"].get_A_inverse_as_sqr() # from dials model, integration step
  permute = sqr((0,0,1,0,1,0,-1,0,0))
  sim_compatible = direct_A*permute # permute columns when post multiplying
  P = crystal_orientation.crystal_orientation(
      sim_compatible,crystal_orientation.basis_type.direct)
  P.show(legend="dials_integrated, C2")
  PR = P.change_basis(CB_OP_C_P)
  PR.show(legend="dials_integrated, primitive setting")
  PRC2 = PR.change_basis(CB_OP_C_P.inverse()) # dials integrated, C2
  cb_op_align = PR.best_similarity_transformation(C,200,1)
  align_PR = PR.change_basis(sqr(cb_op_align))
  align_PR.show(legend="dials_integrated, P1, aligned")
  print("alignment matrix", cb_op_align)
  metric_val = metric(align_PR,C)
  print("Key %d aligned angular offset is %12.9f deg."%(idx, metric_val))
  print("Key %d alignC2 angular offset is %12.9f deg."%(idx, metric(align_PR.change_basis(CB_OP_C_P.inverse()),C2)))
  minimum = metric_val
  ixm = 0
  iym = 0
  best_Ori_C2 = None
  for ix in range(-10,11,2):
    xrotated_Ori = PRC2.rotate_thru((0,0,1),ix*0.01*math.pi/180.)
    for iy in range(-10,11,2):
      yrotated_Ori = xrotated_Ori.rotate_thru((0,1,0),iy*0.01*math.pi/180.)
      new_aligned_Ori = yrotated_Ori.change_basis(CB_OP_C_P).change_basis(sqr(cb_op_align)).change_basis(CB_OP_C_P.inverse())
      grid_metric_val = metric(new_aligned_Ori,C2)
      if grid_metric_val<minimum:
        ixm = ix
        iym = iy
        best_Ori_C2 = new_aligned_Ori
      #print("ix %4d"%ix,"iy %4d"%iy,"grid search angular offset is %12.9f deg."%(grid_metric_val))
      minimum = min(minimum, grid_metric_val)

  print("Key %d minimC2 angular offset is %12.9f deg."%(idx,minimum),"with ix=%d iy=%d"%(ixm,iym))
  best_Ori = best_Ori_C2.change_basis(CB_OP_C_P)
  print("Key %d minimum angular offset is %12.9f deg."%(idx,metric(best_Ori,C)),"with ix=%d iy=%d"%(ixm,iym))
  best_Ori.show(legend="superpower, aligned")
  return dict(superpower_postrefine=best_Ori.direct_matrix(),
              dials_refine=align_PR.direct_matrix(),
              coarse_ground_truth=C.direct_matrix())

def run_sim2smv(ROI,prefix,crystal,spectra,rotation,rank,tophat_spectrum=True,quick=False):
  smv_fileout = prefix + ".img"

  direct_algo_res_limit = 1.7

  wavlen, flux, wavelength_A = next(spectra) # list of lambdas, list of fluxes, average wavelength

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

  #from matplotlib import pyplot as plt
  #plt.plot(flux,"r-")
  #plt.title(smv_fileout)
  #plt.show()

  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=data().get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  CB_OP_C_P = GF.xray_structure.change_of_basis_op_to_primitive_setting() # from C to P, used for ersatz model
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
  #import sys
  #if len(sys.argv)>2:
  #  SIM.seed = -int(sys.argv[2])
  #  print("GOTHERE seed=",SIM.seed)
  #if len(sys.argv)>1:
  #  if sys.argv[1]=="random" : SIM.randomize_orientation()
  SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
  SIM.mosaic_domains = 50  # 77 seconds.  With 100 energy points, 7700 seconds (2 hours) per image
                           # 3000000 images would be 100000 hours on a 60-core machine (dials), or 11.4 years
                           # using 2 nodes, 5.7 years.  Do this at SLAC? NERSC? combination of all?
                           # SLAC downtimes: Tues Dec 5 (24 hrs), Mon Dec 11 (72 hrs), Mon Dec 18 light use, 24 days
                           # mosaic_domains setter must come after mosaic_spread_deg setter
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

  #####  ersatz inserted code (3 lines) for modeling profiles from refined geometry ###
  key = int(prefix.split("_")[-1])
  different_models = superpower_postrefine(key,CB_OP_C_P,old_Amat=Amatrix_rot)
  for key in different_models:
    print (key, sqr(different_models[key]).elems)
  new_Amat = sqr(different_models[model_mode])
  Amatrix_rot = new_Amat # insert the chosen-model Amatrix
  print ("&&& chosen",new_Amat.elems)

  SIM.Amatrix_RUB = Amatrix_rot
  #workaround for failing init_cell, use custom written Amatrix setter
  print("unit_cell_Adeg=",SIM.unit_cell_Adeg)
  print("unit_cell_tuple=",SIM.unit_cell_tuple)
  Amat = sqr(SIM.Amatrix).transpose() # recovered Amatrix from SIM
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
  make_response_plot = response_plot(False,title=prefix)

  for x in range(0,100,2): #len(flux)):
    if flux[x]==0.0:continue
    print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
    CH = channel_pixels(ROI,wavlen[x],flux[x],N,UMAT_nm,Amatrix_rot,GF,output)
    incremental_signal = CH.raw_pixels * crystal.domains_per_crystal
    make_response_plot.append_channel(x,ROI,incremental_signal)
    if x in [26,40,54,68]: # subsample 7096, 7110, 7124, 7138 eV
      print ("----------------------> subsample", x)
      make_response_plot.incr_subsample(x,ROI,incremental_signal)
    make_response_plot.increment(x,ROI,incremental_signal)
    SIM.raw_pixels += incremental_signal;
    CH.free_all()

  message = output.getvalue().split()
  miller = (int(message[4]),int(message[5]),int(message[6]))
  intensity = float(message[9]);

  #SIM.to_smv_format(fileout=prefix + "_intimage_001.img")
  pixels = SIM.raw_pixels
  roi_pixels = pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
  print("Reducing full shape of",pixels.focus(),"to ROI of",roi_pixels.focus())
  SIM.free_all()
  make_response_plot.plot(roi_pixels,miller)
  return dict(roi_pixels=roi_pixels,miller=miller,intensity=intensity,
              channels=make_response_plot.channels)

class response_plot:
  def __init__(self,enable=True,title=""):
    self.enable=enable
    self.signal = None
    self.subsample = None # sample just a few wavelengths
    self.denominator = None # to normalize the sum
    self.title=title
    self.channels = {}
  def append_channel(self, channel, ROI, incr):
    roi_incr = incr[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
    self.channels[channel] = roi_incr
  def increment(self, channel, ROI, incr):
    if not self.enable: return
    roi_incr = incr[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
    signal = channel * roi_incr
    if self.signal is None:
      self.signal = signal;
    else:
      self.signal += signal;
  def incr_subsample(self, channel, ROI, incr):
    if not self.enable: return
    roi_incr = incr[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
    subsample = channel * roi_incr
    if self.subsample is None:
      self.subsample = subsample; self.denominator = roi_incr
    else:
      self.subsample += subsample; self.denominator += roi_incr
  def plot(self,pixels,miller):
        if not self.enable: return
        import numpy as np
        from matplotlib import pyplot as plt

        sig2 = self.signal < 1.
        sig3 = pixels.set_selected(sig2,1.E-200)
        energy_each_pxl_avg = (self.signal/sig3).set_selected(sig2,255.)

        values_1 = pixels # sb.data # ADU
        v0_1 = values_1.set_selected(values_1<=0, 0.)
        v1_1 = v0_1.set_selected(v0_1>255,255)
        v2_1 = (256.-v1_1)/256.
        np_v2_1 = np.ndarray(shape=pixels.focus(), dtype=np.float64, buffer=v2_1.as_numpy_array())

        # this two-d array encodes the color, (0=7070eV, 100=7170eV), now construct 3D RGB array.
        from matplotlib.colors import hsv_to_rgb
        HSV = np.ndarray(shape=(pixels.focus()[0],pixels.focus()[1],3), dtype=np.float64)
        # hue comes directly from the energy:
        HSV[:,:,0] = (energy_each_pxl_avg.as_numpy_array()-20.)/100. #-20 reddens it a bit
                     # red = 0, 7090 eV    blue=0.6, 7150 eV
        HSV[:,:,1] = 1.-np_v2_1
        HSV[:,:,2] = 1.0
        RGB = hsv_to_rgb(HSV)

        sig3ss = self.denominator.deep_copy().set_selected(sig2,1.E-200) # deep copy avoids modifying the pixels array itself
        energyss_each_pxl_avg = (self.subsample/sig3ss).set_selected(sig2,255.)

        values_1ss = self.denominator/pixels
        v0_1ss = values_1ss.set_selected(pixels<=1., 0.)
        np_v2_1ss = np.ndarray(shape=pixels.focus(), dtype=np.float64, buffer=v0_1ss.as_numpy_array())

        HSVss = np.ndarray(shape=(pixels.focus()[0],pixels.focus()[1],3), dtype=np.float64)
        # hue comes directly from the energy:
        HSVss[:,:,0] = (energyss_each_pxl_avg.as_numpy_array()-20.)/100. #-20 reddens it a bit
                     # red = 0, 7090 eV    blue=0.6, 7150 eV
        HSVss[:,:,1] = np_v2_1ss # (self.denominator / pixels).as_numpy_array()
        HSVss[:,:,2] = 1.0
        RGBss = hsv_to_rgb(HSVss)

        ssplot=True # XXXYYY
        if ssplot:
          ax=ax1=ax2=None
          fig = plt.figure(figsize=(8,9))
          fig.suptitle(t="Image %s %s"%(self.title,str(miller)))
          ax = plt.subplot2grid ((3,1),(0,0))
          ax.imshow(np_v2_1, cmap=plt.cm.gray, interpolation="nearest")
          ax1 = plt.subplot2grid ((3,1),(1,0))
          ax1imsh = ax1.imshow(RGB,  interpolation="nearest")
          color_vals = [1.0,0.8,0.6,0.4]
          ax1cb = plt.colorbar(ax1imsh,ticks = color_vals, label = "Energy (eV)")
               # 1.0 = 7090 eV, 0.8 = 7110 eV, 0.6 = 7130 eV, 0.4 = 7150 eV
          ax1cb.ax.set_yticklabels(["%4.0f"%(i) for i in [7090., 7110., 7130., 7150.]])
          ax1cb.ax.set_ylim(1.01,0.39)

          ax2 = plt.subplot2grid ((3,1),(2,0))
          ax2imsh = ax2.imshow(RGBss,  interpolation="nearest")
          color_vals = [1.0,0.8,0.6,0.4]
          ax2cb = plt.colorbar(ax2imsh,ticks = color_vals, label = "Energy (eV)")
               # 1.0 = 7090 eV, 0.8 = 7110 eV, 0.6 = 7130 eV, 0.4 = 7150 eV
          plt.show()

def get_partiality_response(key,one_index,spectra_simulation,ROI,rand_ori,tophat_spectrum=True):

  spectra = spectra_simulation
  crystal = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals

  iterator = spectra.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)

  file_prefix = "key_slow_nonoise_%06d"%key

  pixels = run_sim2smv(ROI,prefix = file_prefix,crystal = crystal,spectra=iterator,rotation=rand_ori,tophat_spectrum=tophat_spectrum,quick=False,rank=0)
  return pixels

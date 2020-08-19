from __future__ import print_function, division
from six.moves import range
from scitbx.array_family import flex
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step5_pad import data
from cctbx import crystal_orientation
from scitbx.matrix import sqr
from simtbx.nanoBragg import shapetype
from simtbx.nanoBragg import nanoBragg
from six.moves import StringIO
import scitbx
import math
from scitbx.matrix import col,sqr
from LS49.work2_for_aca_lsq.util_partiality import channel_pixels
import os
from six.moves import cPickle as pickle

json_glob = os.environ["JSON_GLOB"]
pickle_glob = os.environ["PICKLE_GLOB"]

from LS49 import legacy_random orientations as ersatz_all_orientations

def get_item_from_key_and_glob(key,abc_glob):
      abc_file = abc_glob%key
      with open(abc_file,"rb") as F:
        T = pickle.load(F)
      return T

class differential_roi_manager(object):
  def __init__(self,key,spotlist,spectrum,crystal,allspectra):
    self.allspectra = allspectra # remove later
    from dxtbx.model.experiment_list import ExperimentListFactory
    from six.moves import cPickle as pickle
    E = ExperimentListFactory.from_json_file(json_glob%key,check_format=False)[0] # the dials experiment
    C = E.crystal
    C.show()
    self.data = pickle.load(open(pickle_glob%key,"rb")) # the dials reflections file
    self.gen_fmodel_adapt() # generate Fmodel to obtain CB_OP_C_P
    self.models4 = self.get_idx_rotation_models(key) # alignment between dials refine and coarse ground truth
    metric_P1, metric_C2 = self.get_current_angular_offsets([0.,0.,0.])
    print("""InitialVal Image %06d on %d Bragg spots angular offsets in P1 and C2 (degrees): %8.5f %8.5f"""%(
        key, len(spotlist), metric_P1, metric_C2))
    self.perform_simulations(spectrum,crystal,tophat_spectrum=False)
    # the following is only for debugging:
    #self.model_rotations_and_spots(key,spotlist) # get shoeboxes and differential rotations

  def __del__(self):
    self.SIM.free_all()

  def perform_simulations(self,spectrum, crystal, tophat_spectrum=True):
    #borrow code from util_partiality.
    #def run_sim2smv(ROI,prefix,crystal,spectra,rotation,rank,tophat_spectrum=True,quick=False):

    direct_algo_res_limit = 1.7

    self.wavlen, self.flux, self.wavelength_A = next(spectrum) # list of lambdas, list of fluxes, average wavelength

    if tophat_spectrum:
      sum_flux = flex.sum(self.flux)
      ave_flux = sum_flux/60. # 60 energy channels
      for ix in range(len(self.wavlen)):
        energy = 12398.425 / self.wavlen[ix]
        if energy>=7090 and energy <=7150:
          self.flux[ix]=ave_flux
        else:
          self.flux[ix]=0.

    # use crystal structure to initialize Fhkl array
    self.sfall_main.show_summary(prefix = "Amplitudes used ")
    N = crystal.number_of_cells(self.sfall_main.unit_cell())
    self.crystal = crystal # delete later
    self.N = N # delete later
    SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(N,N,N),
      # workaround for problem with wavelength array, specify it separately in constructor.
      wavelength_A=self.wavelength_A,verbose=0)
    self.SIM = SIM
    SIM.adc_offset_adu = 10 # Do not offset by 40
    SIM.mosaic_spread_deg = 0.05 # interpreted by UMAT_nm as a half-width stddev
    SIM.mosaic_domains = 50  # mosaic_domains setter must come after mosaic_spread_deg setter
#manuscript says 200, June 15 abc_cov was done with 50
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
    self.UMAT_nm = UMAT_nm # delete later

    # get same noise each time this test is run
    SIM.seed = 1
    SIM.oversample=1
    SIM.wavelength_A = self.wavelength_A
    SIM.polarization=1
    # this will become F000, marking the beam center
    SIM.default_F=0
    #SIM.Fhkl=self.sfall_main # no it turns out we don't use these calculations for abc_coverage

    # use local file with (open(something,"wb")) as F:
    with (open("confirm_sfall_P1_7122_amplitudes.pickle","rb")) as F:
      sfall_P1_7122_amplitudes = pickle.load(F)
    SIM.Fhkl = sfall_P1_7122_amplitudes


    SIM.xtal_shape=shapetype.Gauss # both crystal & RLP are Gaussian
    SIM.progress_meter=False
    SIM.show_params()
    # flux is always in photons/s
    SIM.flux=1e12
    SIM.exposure_s=1.0 # so total fluence is e12
    # assumes round beam
    SIM.beamsize_mm=0.003 #cannot make this 3 microns; spots are too intense
    temp=SIM.Ncells_abc
    print("Ncells_abc=",SIM.Ncells_abc)
    SIM.Ncells_abc=temp

  def perform_one_simulation(self):
    ROI = self.ROI
    Amatrix_rot = self.models4["Amat"]
    self.SIM.Amatrix_RUB = Amatrix_rot
    #workaround for failing init_cell, use custom written Amatrix setter
    Amatrecover = sqr(self.SIM.Amatrix).transpose() # recovered Amatrix from SIM
    Ori = crystal_orientation.crystal_orientation(Amatrecover, crystal_orientation.basis_type.reciprocal)
    self.SIM.raw_pixels.fill(0.0) # effectively initializes the data pixels for a new simulation

    self.SIM.seed = 1
    # simulated crystal is only 125 unit cells (25 nm wide)
    # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
    output = StringIO() # open("myfile","w")
    #  make_response_plot = response_plot(False,title=prefix)

    for x in range(0,100,2): #len(flux)):
      if self.flux[x]==0.0:continue
      #print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
      CH = channel_pixels(ROI,self.wavlen[x],self.flux[x],self.N,self.UMAT_nm,Amatrix_rot,self.GF,output)
      incremental_signal = CH.raw_pixels * self.crystal.domains_per_crystal
      #  make_response_plot.append_channel(x,ROI,incremental_signal)
      # if x in [26,40,54,68]: # subsample 7096, 7110, 7124, 7138 eV
      #   print ("----------------------> subsample", x)
      #   make_response_plot.incr_subsample(x,ROI,incremental_signal)
      #make_response_plot.increment(x,ROI,incremental_signal)
      self.SIM.raw_pixels += incremental_signal;
      CH.free_all()

    #message = output.getvalue().split()
    #miller = (int(message[4]),int(message[5]),int(message[6]))
    #intensity = float(message[9]);

    pixels = self.SIM.raw_pixels
    roi_pixels = pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
    #  make_response_plot.plot(roi_pixels,miller)
    #  return dict(roi_pixels=roi_pixels,miller=miller,intensity=intensity,
             # channels=make_response_plot.channels)
    return roi_pixels

  def perform_one_simulation_optimized(self,model,ROI,models4):
    # models4 is a dict of direct_matrices with keys "Amat","Amat_dx", "Amat_dy", "Amat_dz"
    Amatrix_rot = models4[model]
    self.SIM.Amatrix_RUB = Amatrix_rot
    #workaround for failing init_cell, use custom written Amatrix setter
    Amatrecover = sqr(self.SIM.Amatrix).transpose() # recovered Amatrix from SIM
    Ori = crystal_orientation.crystal_orientation(Amatrecover, crystal_orientation.basis_type.reciprocal)
    self.SIM.raw_pixels.fill(0.0) # effectively initializes the data pixels for a new simulation
    temp_fill_zeroes_for_roi = self.SIM.raw_pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]
    channel_summation_roi_only = self.SIM.raw_pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]]

    self.SIM.seed = 1
    # simulated crystal is only 125 unit cells (25 nm wide)
    # amplify spot signal to simulate physical crystal of 4000x larger: 100 um (64e9 x the volume)
    output = StringIO() # open("myfile","w")
    prefix = "key_slow_nonoise"
    from LS49.work2_for_aca_lsq.util_partiality import response_plot
    make_response_plot = response_plot(enable=False,title=prefix)

    self.SIM.region_of_interest = ROI

    from boost_adaptbx.boost.python import streambuf
    self.SIM.printout=True # only for the purpose of getting the P1 Miller index and structure factor
    fast = ROI[1][0] + (ROI[1][1]-ROI[1][0])//2
    slow = ROI[0][0] + (ROI[0][1]-ROI[0][0])//2
    self.SIM.printout_pixel_fastslow=(slow,fast)

    for x in range(0,100,2): #len(flux)):
      if self.flux[x]==0.0:continue
      # initialize ROI-only to zero
      self.SIM.raw_pixels.matrix_paste_block_in_place(temp_fill_zeroes_for_roi,ROI[1][0],ROI[0][0])
      self.SIM.wavelength_A = self.wavlen[x]
      self.SIM.flux = self.flux[x]
      self.SIM.add_nanoBragg_spots_nks(streambuf(output))
      self.SIM.printout = False # only do the printout once through
      incremental_signal = self.SIM.raw_pixels[ROI[1][0]:ROI[1][1], ROI[0][0]:ROI[0][1]] * self.crystal.domains_per_crystal
      make_response_plot.append_channel_full_region(x,incremental_signal)
      if x in [26,40,54,68]: # subsample 7096, 7110, 7124, 7138 eV
        #print ("----------------------> subsample", x)
        make_response_plot.incr_subsample_full_region(x,incremental_signal)
      make_response_plot.increment_full_region(x,incremental_signal)
      channel_summation_roi_only += incremental_signal;

    message = output.getvalue().split()
    miller = (int(message[4]),int(message[5]),int(message[6])) # nanoBragg P1 index
    intensity = float(message[9]);

    make_response_plot.plot(channel_summation_roi_only,miller)
    return dict(roi_pixels=channel_summation_roi_only,miller=miller,intensity=intensity,
                channels=make_response_plot.channels)

  def XXXperform_one_simulation(self,key):
    # method uses the old code as used in abc_coverage
    #def get_partiality_response(key,one_index,spectra_simulation,ROI,rand_ori,tophat_spectrum=True):

    spectra = self.allspectra
    from LS49.sim.step5_pad import microcrystal
    crystal = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0) # assume smaller than 10 um crystals

    iterator = spectra.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)

    file_prefix = "key_slow_nonoise_%06d"%key
    Amatrix_rot = self.models4["Amat"]
    from LS49.work2_for_aca_lsq.util_partiality import run_sim2smv
    pixels = run_sim2smv(self.ROI,prefix = file_prefix,crystal = crystal,spectra=iterator,rotation=Amatrix_rot,tophat_spectrum=False,quick=False,rank=0)
    return pixels["roi_pixels"]


  def model_rotations_and_spots(self,key,spotlist):
    #ersatz import of other refinements
    abc_dials_refine = get_item_from_key_and_glob(key,os.environ["ABC_GLOB_A"])
    abc_coarse_truth = get_item_from_key_and_glob(key,os.environ["ABC_GLOB_C"])

    M = self.data["miller_index"]
    for model in self.models4:
      if model != "Amat": continue
      print("\n",model)
      for ispot, spot in enumerate(spotlist):
        ('a', 'asu_idx_C2_setting', 'bkgrd_a', 'channels', 'compute_functional_and_gradients', 'image_no', 'n', 'orig_idx_C2_setting',
         'print_step', 'roi', 'sb_data', 'simtbx_P1_miller', 'simtbx_intensity_7122', 'x')
        S = (spot.orig_idx_C2_setting) # happens to be the coarse ground truth spot
        abc_dials_refine_spot = [spotL for spotL in abc_dials_refine if spotL.orig_idx_C2_setting==S][0]
        abc_coarse_truth_spot = [spotL for spotL in abc_coarse_truth if spotL.orig_idx_C2_setting==S][0]

        idx = M.first_index(S)
        shoe = self.data["shoebox"][idx]
        B = shoe.bbox
        ROI = ((B[0],B[1]),(B[2],B[3]))
        print ("C2",S,ROI)

        from LS49.ML_push.shoebox_troubleshoot import pprint3,pprint
        #pprint3 (shoe.data) # shoe.data is identical to spot.sb_data, in case you were wondering
#        pprint (abc_coarse_truth_spot.roi)
#        pprint (abc_dials_refine_spot.roi)
        #new_calc = self.perform_one_simulation()
        #pprint (new_calc)
        new_calc2_dict = self.perform_one_simulation_optimized(model,ROI,self.models4)
        pprint (new_calc2_dict["roi_pixels"])

  def gen_fmodel_adapt(self):
    direct_algo_res_limit = 1.7
    self.GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=data(
         ).get("pdb_lines"),algorithm="fft",wavelength=7122)
    self.CB_OP_C_P = self.GF.xray_structure.change_of_basis_op_to_primitive_setting() # from C to P, for ersatz model
    self.GF.set_k_sol(0.435)
    self.GF.make_P1_primitive()
    self.sfall_main = self.GF.get_amplitudes()

  def get_idx_rotation_models(self,idx):
    rotation = sqr(ersatz_all_orientations()[idx])
    Amat = (rotation * sqr(self.sfall_main.unit_cell().orthogonalization_matrix())).transpose()
    from LS49.work_pre_experiment.post5_ang_misset import get_item
    from LS49.ML_push.exploratory_missetting import metric
    R = get_item(idx)
    print ("coarse ground truth with index",idx)
    C = crystal_orientation.crystal_orientation(
      Amat,crystal_orientation.basis_type.direct)
    C.show(legend="ground truth, P1")
    C2 = C.change_basis(self.CB_OP_C_P.inverse())
    C2.show(legend="ground truth, C2")
    direct_A = R["integrated_crystal_model"].get_A_inverse_as_sqr() # from dials model, integration step
    permute = sqr((0,0,1,0,1,0,-1,0,0))
    sim_compatible = direct_A*permute # permute columns when post multiplying
    P = crystal_orientation.crystal_orientation(
      sim_compatible,crystal_orientation.basis_type.direct)
    P.show(legend="dials_integrated, C2")
    PR = P.change_basis(self.CB_OP_C_P)
    PR.show(legend="dials_integrated, primitive setting")
    PRC2 = PR.change_basis(self.CB_OP_C_P.inverse()) # dials integrated, C2
    cb_op_align = PR.best_similarity_transformation(C,200,1)
    align_PR = PR.change_basis(sqr(cb_op_align))
    align_PR.show(legend="dials_integrated, P1, aligned")
    print("alignment matrix", cb_op_align)
    metric_val = metric(align_PR,C)
    print("Key %d aligned angular offset is %12.9f deg."%(idx, metric_val))
    print("Key %d alignC2 angular offset is %12.9f deg."%(idx, metric(align_PR.change_basis(self.CB_OP_C_P.inverse()),C2)))
    # coarse, dials crystal orientation models = C, align_PR
    # apply Rotx:
    align_PR_dx = align_PR.rotate_thru((1.0,0.0,0.0), math.pi* 0.01/180.)
    align_PR_dy = align_PR.rotate_thru((0.0,1.0,0.0), math.pi* 0.01/180.)
    align_PR_dz = align_PR.rotate_thru((0.0,0.0,1.0), math.pi* 0.01/180.)
    self.dials_integrated_P1_aligned = align_PR; self.ground_truth_P1 = C; self.ground_truth_C2 = C2
    return (dict(Amat=align_PR.direct_matrix(),Amat_dx=align_PR_dx.direct_matrix(),
            Amat_dy=align_PR_dy.direct_matrix(), Amat_dz=align_PR_dz.direct_matrix()))

  def get_incremented_rotation_models(self,rotxyz): # pass in the x,y,z rotations in units of 0.01 degree
    from LS49.ML_push.exploratory_missetting import metric

    updated_mat = self.dials_integrated_P1_aligned.rotate_thru((1.0,0.0,0.0), rotxyz[0] * math.pi* 0.01/180.
                                                 ).rotate_thru((0.0,1.0,0.0), rotxyz[1] * math.pi* 0.01/180.
                                                 ).rotate_thru((0.0,0.0,1.0), rotxyz[2] * math.pi* 0.01/180.)

    metric_val = metric(updated_mat,self.ground_truth_P1)
    print("minimization step, aligned angular offset is %12.9f deg."%(metric_val))
    new_dx = updated_mat.rotate_thru((1.0,0.0,0.0), math.pi* 0.01/180.)
    new_dy = updated_mat.rotate_thru((0.0,1.0,0.0), math.pi* 0.01/180.)
    new_dz = updated_mat.rotate_thru((0.0,0.0,1.0), math.pi* 0.01/180.)
    return (dict(Amat=updated_mat.direct_matrix(),Amat_dx=new_dx.direct_matrix(),
            Amat_dy=new_dy.direct_matrix(), Amat_dz=new_dz.direct_matrix()))

  def get_current_angular_offsets(self,rotxyz): # pass in the x,y,z rotations in units of 0.01 degree
    from LS49.ML_push.exploratory_missetting import metric

    updated_mat = self.dials_integrated_P1_aligned.rotate_thru((1.0,0.0,0.0), rotxyz[0] * math.pi* 0.01/180.
                                                 ).rotate_thru((0.0,1.0,0.0), rotxyz[1] * math.pi* 0.01/180.
                                                 ).rotate_thru((0.0,0.0,1.0), rotxyz[2] * math.pi* 0.01/180.)

    metric_P1 = metric(updated_mat,self.ground_truth_P1)
    metric_C2 = metric(updated_mat.change_basis(self.CB_OP_C_P.inverse()) ,self.ground_truth_C2)

    return metric_P1,metric_C2

"""
work out the basic code for
OK 1) getting CB_OP
OK 2) Getting sfall_main
OK 3) getting the coarse ground truth
OK 4) getting the dials refine
OK 5) applying rotational perturbations to dials refine
6) performing 7150 eV simulations for all three orientations (GT, RXGT, RYGT)
"""

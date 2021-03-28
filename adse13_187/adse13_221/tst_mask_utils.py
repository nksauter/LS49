from __future__ import division
from dials.array_family import flex
import pickle
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.development.timers import Profiler
import math

def top_75_iterator():
  import re
  from LS49.adse13_187.case_data import lookup_repo
  for ikey, key in enumerate(lookup_repo):
    assert ikey==key
    item = lookup_repo[key]
    match = re.search("shot([0-9]*)",item)
    event_idx = int(match.groups()[0])
    yield event_idx

class mask_manager:
  def __init__(self,trusted_mask,refl_table,expt):
    self.trusted_mask = trusted_mask
    self.refl_table = refl_table
    self.expt = expt

  def mask_introspection(self):
    for midx in range(256):
      print(midx, (self.trusted_mask[midx].count(True)), (self.trusted_mask[midx].count(False)))
    #conclusion:  mask True==good pixels; mask False==bad pixels
    # invert mask:  True==bad; False==good
    # on the invert mask, we can set the shoebox pixels to False

  def get_trusted_and_refl_mask(self):
    mask = self.trusted_mask
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)

    self.shoebox_mask=[]
    for midx in range(256):
      self.shoebox_mask.append(flex.bool(flex.grid(mask[midx].focus()), False))
      # should be a False field

    N_true = 0
    N_bad_pizel = 0
    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      for islow in range(max(0,bbox[2]-3), min(slow_size,bbox[3]+3)):
        for ifast in range(max(0,bbox[0]-3), min(fast_size,bbox[1]+3)):
          false_field[islow*slow_size + ifast]=True
          N_true += 1
          if mask[ipanel][islow*slow_size + ifast] is False: N_bad_pizel+=1
    self.resultant = []
    for midx in range(256):
      newmask = mask[midx].__and__(self.shoebox_mask[midx])
      self.resultant.append(newmask)
      if False and midx<4: print ("Panel",midx, "good px", mask[midx].count(True), "shoebox px",
             self.shoebox_mask[midx].count(True), "result", newmask.count(True))
    print (N_true,"pixels were visited in the %d shoeboxes (with borders)"%size)
    print (N_bad_pizel,"of these were bad pixels")

  def refl_analysis(self,dials_model):
    """This function sets up some data structures (spots_*) allowing us to index into the spots
    and pixels of interest.  These will be repeatedly used during parameter refinement
    to calculate target function and intermediate statistics.
    """
    Z = self.refl_table
    indices = Z['miller_index']
    expts = ExperimentListFactory.from_json_file(dials_model,
                                              check_format=False)
    self.dials_model=expts[0]
    CRYS = self.dials_model.crystal
    UC = CRYS.get_unit_cell()
    strong_resolutions = UC.d(indices)
    order = flex.sort_permutation(strong_resolutions, reverse=True)
    Z["spots_order"] = order
    self.spots_pixels = flex.size_t()
    spots_offset = flex.int(len(order),-1)
    spots_size = flex.int(len(order),-1)

    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    N_visited = 0; N_bad = 0
    for oidx in range(len(order)): #loop through the shoeboxes in correct order
      sidx = order[oidx] # index into the Miller indices
      ipanel = P[sidx]
      slow_size = 254
      fast_size = 254
      panel_size=slow_size*fast_size
      bbox = S[sidx].bbox
      first_position = spots_offset[sidx] = self.spots_pixels.size()
      for islow in range(max(0,bbox[2]-3), min(slow_size,bbox[3]+3)):
        for ifast in range(max(0,bbox[0]-3), min(fast_size,bbox[1]+3)):
          value = self.trusted_mask[ipanel][islow*slow_size + ifast]
          N_visited += 1
          if value: self.spots_pixels.append(ipanel*panel_size+islow*slow_size+ifast)
          else: N_bad+=1
      spot_size = spots_size[sidx] = self.spots_pixels.size() - first_position
    Z["spots_offset"] = spots_offset
    Z["spots_size"] = spots_size
    print (N_visited,"pixels were visited in the %d shoeboxes (with borders)"%len(order))
    print (N_bad,"of these were bad pixels, leaving %d in target"%(len(self.spots_pixels)))

  def simple_rmsd(self,calc_data="xyzcal.px",plot=False):
    """Function does a rudimentary plot of model vs. experimental spot position, and optionally
    a plot of deviation vs. Bragg angle.
    """
    Z = self.refl_table
    devs = flex.double()
    sqdevs = flex.double()
    for oidx in range(len(Z["spots_order"])): #loop through the shoeboxes in correct order
      sidx = Z["spots_order"][oidx] # index into the Miller indices
      calc = Z[calc_data][sidx][0:2]
      obs  = Z['xyzobs.px.value'][sidx][0:2]
      sqdev= (calc[0]-obs[0])**2 + (calc[1]-obs[1])**2
      dev  = math.sqrt(sqdev)
      devs.append(dev)
      sqdevs.append(sqdev)
      #print("%20s %6.2fpx"%(Z["miller_index"][sidx], dev))
    print ("The rmsd is %6.2f px"%math.sqrt(flex.mean(sqdevs)))
    if plot:
      from matplotlib import pyplot as plt
      plt.plot(range(len(devs)),devs)
      running_range = range(15,len(devs),15)
      plt.plot(running_range, [flex.mean(devs[I-15:I]) for I in running_range], "r-",
         label="RMSD=%6.2fpx"%math.sqrt(flex.mean(sqdevs)))
      plt.title("Model vs. Experimental Bragg spot position")
      plt.xlabel("Spots ordered by increasing Bragg angle →")
      plt.ylabel("Deviation in pixels")
      plt.legend(loc='upper left')
      plt.show()

  def plot_pixel_histograms(self):
    exp_data = self.expt.imageset.get_raw_data(0) # experimental data
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)
    all_pixels = flex.double()
    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      islow_limits = (max(0,bbox[2]-3), min(slow_size,bbox[3]+3))
      ifast_limits = (max(0,bbox[0]-3), min(fast_size,bbox[1]+3))
      for islow in range(islow_limits[0], islow_limits[1]):
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          all_pixels.append(exp_data[ipanel][islow*slow_size + ifast])
    from matplotlib import pyplot as plt
    #plt.hist(all_pixels, len(all_pixels)//100, range=(-10,20000), facecolor="orange")
    plt.hist(all_pixels, len(all_pixels)//100, facecolor="orange")
    #plt.ylim((-10,500))
    plt.ylim((-10,100))
    plt.show()

  def resultant_mask_to_file(self,file_name):
    all_done = tuple(self.resultant)
    with open(file_name,"wb") as M:
      pickle.dump(all_done,M)

  @classmethod
  def from_files(cls, trusted_mask_file, refl_file, expt_file):
    with open(trusted_mask_file,"rb") as M:
      mask = pickle.load(M)
    refl_table = flex.reflection_table.from_file(refl_file)
    expts = ExperimentListFactory.from_json_file(expt_file,
                                              check_format=True)
    return cls(mask,refl_table,expts[0])

  def get_lunus_repl(self):
    P = Profiler("LUNUS")
    # first get the lunus image
    from lunus.command_line.filter_peaks import get_image_params
    from lunus import LunusDIFFIMAGE

    imageset = self.expt.imageset
    data = imageset[0]
    assert isinstance(data, tuple) # assume a tuple of flex::double over detector panels
    # Instantiate a LUNUS diffraction image
    A = LunusDIFFIMAGE(len(data))

    # Populate the image with multipanel data

    for pidx in range(len(data)):
      A.set_image(pidx,data[pidx])

    # Define the LUNUS image parameters
    deck = '''
#lunus input deck
#punchim_xmin=1203
#punchim_ymin=1250
#punchim_xmax=2459
#punchim_ymax=1314
#windim_xmin=100
#windim_ymin=100
#windim_xmax=2362
#windim_ymax=2426
#thrshim_min=0
#thrshim_max=50
modeim_bin_size=1
modeim_kernel_width=15
'''
    image_params = get_image_params(imageset)

    # Set the LUNUS image parameters
    for pidx in range(len(image_params)):
        deck_and_extras = deck+image_params[pidx]
        A.LunusSetparamsim(pidx,deck_and_extras)
    A.LunusModeim()

    # Get the processed image
    lunus_filtered_data = flex.double()
    assert len(data)==256 # Jungfrau
    for pidx in range(len(data)):
      aye_panel = A.get_image_double(pidx)
      assert aye_panel.focus() == (254, 254)
      lunus_filtered_data.extend( aye_panel.as_1d() )
    lunus_filtered_data.reshape(flex.grid((256,254,254)))
    self.lunus_filtered_data = lunus_filtered_data.as_numpy_array()

  def get_image_res_data(self):
    if True: # params.write_experimental_data:
      self.exp_data = self.expt.imageset.get_raw_data(0) # why a different access pattern? Are the data different?
      assert len(self.exp_data) == 256 # Jungfrau panels
      self.exp_data = [self.exp_data[pid].as_numpy_array() for pid in range(len(self.exp_data))]

  def write_hdf5(self,filenm):
    # then write the data
    from simtbx.nanoBragg import utils
    if True: # params.write_output:
      img_sh = self.lunus_filtered_data.shape
      assert img_sh == (256,254,254)
      num_output_images = 6 # 1 + int(params.write_experimental_data)
      print("Saving exascale output data of shape", img_sh)
      beam_dict = self.expt.beam.to_dict()
      det_dict = self.expt.detector.to_dict()
      try:
        beam_dict.pop("spectrum_energies")
        beam_dict.pop("spectrum_weights")
      except Exception: pass

      with utils.H5AttributeGeomWriter(filenm,
                                image_shape=img_sh, num_images=num_output_images,
                                detector=det_dict, beam=beam_dict,
                                detector_and_beam_are_dicts=True) as writer:
        #Output 1.  Lunus pixel-assimilated image
        writer.add_image(self.lunus_filtered_data)

        #Output 2.  In-memory modify the Lunus image, with 1st-order Taylor shoeboxes
        self.modify_shoeboxes()
        writer.add_image(self.lunus_filtered_data)

        if True: # params.write_experimental_data:
          sim_mock = self.simulation_mockup(self.exp_data)

          #Output no. ersatz simulation
          nanobragg_sim = self.ersatz_MCMC()
          #writer.add_image(nanobragg_sim) #hook to produce actual simulation, bypass for now

          #Output 3. analyze proposal and add background
          bragg_plus_background = self.reusable_rmsd(proposal=nanobragg_sim, label="ersatz_mcmc")
          writer.add_image(bragg_plus_background)

          #Output 4. renormalize the proposal
          renormalize_bragg_plus_background = self.reusable_rmsd(proposal=self.renormalize(
            proposal=nanobragg_sim,proposal_label="ersatz_mcmc",ref_label="spots_mockup"),
            label="renormalize_mcmc")
          writer.add_image(renormalize_bragg_plus_background)

          #Output 5. Mockup simulation laid on top of 1st-Taylor background
          writer.add_image(sim_mock)

          #Output 6. Experimental res-data
          writer.add_image(self.exp_data)
        print("Saved output to file %s" % (filenm))

  def ersatz_MCMC(self):
    from LS49.adse13_187.adse13_221.mcmc_class import MCMC_manager
    self.MCMC = MCMC_manager()
    self.MCMC.get_amplitudes(self.dials_model, self.refl_table)
    return self.MCMC.job_runner() # returns simulated image as numpy array

  def per_shoebox_whitelist_iterator(self, sidx):
    """given a shoebox id, iterate through all its whitelisted pixels"""
    Z = self.refl_table
    SOFF = Z["spots_offset"]
    SSIZ = Z["spots_size"]
    slow_size = 254
    panel_size = 254 * 254
    for idxpx in self.spots_pixels[SOFF[sidx]:SOFF[sidx]+SSIZ[sidx]]:
      ipanel = idxpx//panel_size; panelpx = idxpx%panel_size
      islow = panelpx//slow_size; ifast = panelpx%slow_size
      yield ipanel, islow, ifast

  def simulation_mockup(self,experimental):
    """Function creates a mockup simulated image consisting of:
      background = self.lunus_filtered_data + Bragg = (experimental - smooth backgrd)
    Provides an alternate view of the background vs. experiment data.  It zeroes out the
    pixels from the shoeboxes of interest,  but then adds the experimental pixels back.
    Here "experimental pixels" means res-data minus planar-fit-lunus-model.
    Function has the more expansive purpose of analyzing the data & lunus background model
    and storing statistics: the data center of mass in the shoebox, and the data sum, to be
    used later for normalizing the model in each shoebox.
    """
    mockup_ctr_of_mass = flex.vec3_double()
    mockup_shoebox_sum = flex.double()
    import copy
    mockup_simulation = copy.deepcopy(self.lunus_filtered_data)
    from scitbx.matrix import col
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      SUM_VEC = col((0.,0.))
      SUM_wt = 0.
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        model_value = ( experimental[ipanel][islow,ifast] -
                        self.lunus_filtered_data[ipanel,islow,ifast])
        SUM_VEC = SUM_VEC + float(model_value) * col((float(islow),float(ifast)))
        SUM_wt += model_value
        mockup_simulation[ipanel,islow,ifast] = experimental[ipanel][islow,ifast]
        #mockup_simulation[ipanel,islow,ifast] = model_value
      c_o_m = SUM_VEC/SUM_wt
      # there is a half pixel offset in our understanding of position
      mockup_ctr_of_mass.append((c_o_m[1]+0.5,c_o_m[0]+0.5,0.0))
      mockup_shoebox_sum.append(SUM_wt)
    self.refl_table["spots_mockup_xyzcal.px"] = mockup_ctr_of_mass
    self.refl_table["spots_mockup_shoebox_sum"] = mockup_shoebox_sum
    self.simple_rmsd(calc_data="spots_mockup_xyzcal.px",plot=False)
    return mockup_simulation

  def reusable_rmsd(self,proposal,label):
    """Function analyzes proposal data consisting of proposed Bragg spots
    Function has the more expansive purpose of analyzing the data
    and storing statistics: the data center of mass in the shoebox, and the data sum, to be
    used later for normalizing the model in each shoebox.
    """
    proposal_ctr_of_mass = flex.vec3_double()
    proposal_shoebox_sum = flex.double()
    import copy
    mockup_simulation = copy.deepcopy(self.lunus_filtered_data)
    from scitbx.matrix import col
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      SUM_VEC = col((0.,0.))
      SUM_wt = 0.
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        proposal_value = proposal[ipanel][islow,ifast]
        SUM_VEC = SUM_VEC + float(proposal_value) * col((float(islow),float(ifast)))
        SUM_wt += proposal_value
        mockup_simulation[ipanel,islow,ifast] += proposal_value
      c_o_m = SUM_VEC/SUM_wt
      # there is a half pixel offset in our understanding of position
      proposal_ctr_of_mass.append((c_o_m[1]+0.5,c_o_m[0]+0.5,0.0))
      proposal_shoebox_sum.append(SUM_wt)
    self.refl_table[label+"_xyzcal.px"] = proposal_ctr_of_mass
    self.refl_table[label+"_shoebox_sum"] = proposal_shoebox_sum
    self.simple_rmsd(calc_data=label+"_xyzcal.px",plot=True)
    return mockup_simulation

  def renormalize(self,proposal,proposal_label,ref_label):
    """Takes one proposal and returns a new one, with each spot adjusted by a
    separate scale factor that represents the ratio of experiment::proposal
    """
    import copy
    renormalized = copy.deepcopy(proposal)
    all_scales = self.refl_table[ref_label+"_shoebox_sum"]/self.refl_table[proposal_label+"_shoebox_sum"]
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      scale_factor = all_scales[sidx]
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        renormalized[ipanel,islow,ifast] *= scale_factor
    return renormalized

  def modify_shoeboxes(self, verbose=False): # and printing the shoeboxes in verbose mode
    exp_data = self.expt.imageset.get_raw_data(0) # experimental data
    Z = self.refl_table
    P = panels = Z['panel']
    S = shoeboxes = Z['shoebox']
    size = len(Z)

    for sidx in range(size): #loop through the shoeboxes
      ipanel = P[sidx]
      false_field = self.shoebox_mask[ipanel]
      slow_size = false_field.focus()[0]
      fast_size = false_field.focus()[1]
      bbox = S[sidx].bbox
      islow_limits = (max(0,bbox[2]-3), min(slow_size,bbox[3]+3))
      ifast_limits = (max(0,bbox[0]-3), min(fast_size,bbox[1]+3))
      if verbose:
       # print out the res-data
       for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = exp_data[ipanel][islow*slow_size + ifast]
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
       print()
       # print out the trusted mask
       flag=True
       for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = float(int(self.resultant[ipanel][islow*slow_size + ifast]))
          if value==False:flag=False
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
       print()
       # print out the lunus-repl shoebox
       for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = self.lunus_filtered_data[ipanel,islow,ifast]
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
       print()
      # now create a 2nd-order fit to the data.  First implementation, no weighting.
      from LS49.adse13_187.adse13_221.smooth_fit import replacement_pixels
      FIT=replacement_pixels(self, ipanel, islow_limits, ifast_limits, shoebox=S[sidx])
      # print out the fit shoebox
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = FIT.model_T(islow,ifast)
          self.lunus_filtered_data[ipanel,islow,ifast]=value # reset lunus array
          if verbose: print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        if verbose: print(" =%6.0f"%(fast_sum/fast_count))
      if verbose: print()
      if verbose:
       # print out the data minus background-fit shoebox
       for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = exp_data[ipanel][islow*slow_size + ifast]-FIT.model_T(islow,ifast)
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
       print("---")
       #input()

def multiple_cases():
  def get_any_case():
    from LS49.adse13_187.adse13_221.ad_hoc_run795_lookup import conversion
    for idx,item in enumerate(top_75_iterator()):
      if idx>0: exit() # quick check the first one
      run_no = 795 # only look at run number 795
      class Empty: pass
      E = Empty()
      E.trusted_mask_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      E.expt_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_%04d.expt"%item
      E.dials_model = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_%04d.expt"%item
      E.out_file = "top_event_%04d.mask"%idx
      E.hdf5_file = "top_exa_%04d.hdf5"%idx
      E.refl_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_%04d.refl"%(
        conversion[run_no][item]
        )
      yield E
  for E in get_any_case():
    print(E.expt_file)
    print(E.refl_file)
    M = mask_manager.from_files(E.trusted_mask_file, E.refl_file, E.expt_file)
    M.get_trusted_and_refl_mask()
    M.refl_analysis(E.dials_model)
    M.simple_rmsd()
    #M.plot_pixel_histograms()
    M.get_lunus_repl() # compute lunus-repl image.
    M.get_image_res_data()
    M.write_hdf5(E.hdf5_file) # write res-data and lunus-repl to an HDF5 file; includes simulation_mockup, for now.
    #M.modify_shoeboxes()
    M.resultant_mask_to_file(E.out_file)

def single_case():
  def get_case_1():
      class Empty: pass
      E = Empty()
      E.trusted_mask_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/trusted_Py3.mask"
      E.refl_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_0309.refl"
      E.expt_file = "/global/cscratch1/sd/nksauter/adse13_187/bernina/split_c/split_0648.expt"
      E.out_file = "XXXevent_648.mask"
      E.hdf5_file = "exa_648.hdf5"
      return E
  E = get_case_1()
  M = mask_manager.from_files(E.trusted_mask_file, E.refl_file, E.expt_file)
  M.get_trusted_and_refl_mask()
  # compute lunus-repl image.  write res-data and lunus-repl to an HDF5 file
  M.get_lunus_repl(E.hdf5_file)
  M.resultant_mask_to_file(E.out_file)

if __name__ == "__main__":
  #single_case()
  multiple_cases()

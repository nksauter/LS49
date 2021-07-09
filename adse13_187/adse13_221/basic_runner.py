from __future__ import division
import os, copy
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from matplotlib import pyplot as plt
import math
from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt

class MCMC_manager:
  def get_amplitudes(self, dials_model, refl_table, test_without_mpi=True):
    D = dials_model
    R = refl_table
    from cctbx.crystal import symmetry
    from cctbx.miller import array,set as miller_set
    uc = D.crystal.get_unit_cell()
    sg = D.crystal.get_space_group()
    MS = miller_set(
         symmetry(unit_cell=uc, space_group=sg), anomalous_flag=True,
         indices=R["miller_index"].select(R["spots_order"]))
    self.amplitudes = array(MS,
         data=flex.sqrt(R["spots_mockup_shoebox_sum"].select(R["spots_order"]))
         )

    from simtbx.gpu import gpu_energy_channels
    self.gpu_channels_singleton = gpu_energy_channels (
        deviceId = 0 ) # determine device by rank id later

  def set_whitelist(self,value):
    self.relevant_whitelist_order = value

from LS49.adse13_187.adse13_221.lunus_wrap import mask_manager
class basic_run_manager(mask_manager):

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
    rmsd = math.sqrt(flex.mean(sqdevs))
    print ("The rmsd is %6.2f px"%rmsd)
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
    return rmsd

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

  def Z_statistics(self, experiment, model, plot=False):
    keV_per_photon = ENERGY_CONV/1000./self.expt.beam.get_wavelength()
    Z_plot = []
    sigma_pixel = []
    import numpy as np
    for x in range(256):
      exp_panel = self.view["exp_data"][x]
      abs_exp_panel = np.abs(exp_panel)
      abs_exp_panel_photons = abs_exp_panel/keV_per_photon # convert to photons
      poisson_noise_sigma = np.sqrt(abs_exp_panel_photons)
      poisson_noise_sigma = np.where(poisson_noise_sigma==0., 1., poisson_noise_sigma)
      sigma_pixel.append(poisson_noise_sigma)
      diff_panel_photons = (model[x] - experiment[x])/keV_per_photon
      offset_Z = (diff_panel_photons/poisson_noise_sigma)*0.1 + 1.0
      Z_plot.append(offset_Z)

    proposal_shoebox_mean_Z = flex.double()
    proposal_shoebox_sigma_Z = flex.double()
    all_Z_values = flex.double()
    sauter_eq_15_likelihood = flex.double()
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      shoebox_Z_values = flex.double()
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        std_dev_denominator_photons = sigma_pixel[ipanel][islow,ifast]
        diff_panel_photons = (model[ipanel][islow,ifast] - experiment[ipanel][islow,ifast])/keV_per_photon
        Z = (diff_panel_photons/std_dev_denominator_photons)
        shoebox_Z_values.append(Z)
        all_Z_values.append(Z)
        sauter_eq_15_likelihood.append( (model[ipanel][islow,ifast]/keV_per_photon) -
          (experiment[ipanel][islow,ifast]/keV_per_photon) * math.log(model[ipanel][islow,ifast]/keV_per_photon))
      stats = flex.mean_and_variance(shoebox_Z_values)
      proposal_shoebox_mean_Z.append( stats.mean() )
      proposal_shoebox_sigma_Z.append( stats.unweighted_sample_standard_deviation() )
    print("proposal negative log likelihood %10f"%(flex.sum(sauter_eq_15_likelihood)))
    stats = flex.mean_and_variance(all_Z_values)
    mnz = stats.mean()
    sgz = stats.unweighted_sample_standard_deviation()
    print("proposal mean Z=%.2f, sigma Z=%.2f"%(mnz, sgz))
    if plot:
      from matplotlib import pyplot as plt
      plt.plot(range(len(self.refl_table)),
               proposal_shoebox_mean_Z.select(self.refl_table["spots_order"]),"r-",
               label="mean Z (all=%.2f)"%(mnz))
      plt.plot(range(len(self.refl_table)),
               proposal_shoebox_sigma_Z.select(self.refl_table["spots_order"]),
               label="std_dev Z (all=%.2f)"%(sgz))
      plt.title("Z_distribution in each shoebox")
      plt.xlabel("Spots ordered by increasing Bragg angle →")
      plt.ylabel("Z-value")
      plt.legend(loc='upper right')
      plt.show()
    return Z_plot

  def ersatz_MCMC(self, variable_params):
    from LS49.adse13_187.adse13_221.case_run import case_job_runner
    class ersatz(MCMC_manager, case_job_runner): pass
    self.MCMC = ersatz()
    self.MCMC.get_amplitudes(self.dials_model, self.refl_table)
    relevant_whitelist_order = flex.size_t()
    for sidx in range(len(self.refl_table["spots_offset"])):#loop through the shoeboxes
      for pidx in range(self.refl_table["spots_offset"][sidx],
                        self.refl_table["spots_offset"][sidx]+self.refl_table["spots_size"][sidx]):
        relevant_whitelist_order.append(self.spots_pixels[pidx])
    self.MCMC.set_whitelist(relevant_whitelist_order)

    modality = "job"
    return self.MCMC.job_runner(expt=self.expt, alt_expt=self.dials_model, params=variable_params,
      mask_array = self.monolithic_mask_whole_detector_as_1D_bool
      ) # returns simulated image as numpy array

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
    mockup_simulation = copy.deepcopy(self.view["plane_shoeboxes"])
    from scitbx.matrix import col
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      SUM_VEC = col((0.,0.))
      SUM_wt = 0.
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        model_value = ( experimental[ipanel][islow,ifast] -
                        self.view["plane_shoeboxes"][ipanel,islow,ifast])
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
    mockup_simulation = copy.deepcopy(self.view["plane_shoeboxes"])
    from scitbx.matrix import col
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      SUM_VEC = col((0.,0.))
      SUM_wt = 0.
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        proposal_value = max(0.1,proposal[ipanel][islow,ifast])
          #workaround for the "Paley" bug. If spot is not predicted, give it some nonzero intensity
        SUM_VEC = SUM_VEC + float(proposal_value) * col((float(islow),float(ifast)))
        SUM_wt += proposal_value
        mockup_simulation[ipanel,islow,ifast] += proposal_value
      c_o_m = SUM_VEC/SUM_wt
      # there is a half pixel offset in our understanding of position
      proposal_ctr_of_mass.append((c_o_m[1]+0.5,c_o_m[0]+0.5,0.0))
      proposal_shoebox_sum.append(SUM_wt)
    self.refl_table[label+"_xyzcal.px"] = proposal_ctr_of_mass
    self.refl_table[label+"_shoebox_sum"] = proposal_shoebox_sum
    self.simple_rmsd(calc_data=label+"_xyzcal.px",plot=False) # toggle for plotting
    return mockup_simulation

  def renormalize(self,proposal,proposal_label,ref_label):
    """Takes one proposal and returns a new one, with each spot adjusted by a
    separate scale factor that represents the ratio of experiment::proposal
    """
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
    self.view["plane_shoeboxes"] = copy.deepcopy(self.view["lunus_filtered_data"])
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
          value = self.view["plane_shoeboxes"][ipanel,islow,ifast]
          print("%6.0f"%value, end="")
          fast_count+=1
          fast_sum+=value
        print(" =%6.0f"%(fast_sum/fast_count))
       print()
      # now create a 2nd-order fit to the data.  First implementation, no weighting.
      from LS49.adse13_187.adse13_221.smooth_fit import replacement_pixels
      FIT=replacement_pixels(self.view["lunus_filtered_data"],
                             ipanel, islow_limits, ifast_limits, shoebox=S[sidx])
      # print out the fit shoebox
      for islow in range(islow_limits[0], islow_limits[1]):
        fast_count=0
        fast_sum=0
        for ifast in range(ifast_limits[0], ifast_limits[1]):
          value = FIT.model_T(islow,ifast)
          self.view["plane_shoeboxes"][ipanel,islow,ifast]=value # reset lunus array
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

def generate_phil_scope():
  from iotbx.phil import parse
  master_phil="""
    include scope LS49.adse13_187.adse13_221.lunus_wrap.phil_scope
    cryst = ""
      .type = path
      .help = The dials refine model, specifically containing an updated dxtbx.crystal model
    model
      .help = Namespace to pass variable model parameters to the nanoBragg simulation
      {
      mosaic_spread {
        value = 0.01
        .type = float(value_min = 0)
        .help = half-width mosaic rotation in degrees, assuming isotropic model
      }
      Nabc {
        refine = True
          .type = bool
        value = (50,50,50)
          .type = ints(size=3, value_min=2)
          .help = domain size along the a,b, and c axes expressed in unit cells
        sigmas = (10,10,10)
          .type = floats(size=3, value_min=0)
        hyperparameter = 0.2
          .type = float(value_min=0.001, value_max=0.9)
          .help = The allowable range of proposal values as a multiplier of current value.
      }
      }
  """
  return parse(master_phil, process_includes=True)
phil_scope = generate_phil_scope()

def parse_input():
  # The script usage
  import libtbx.load_env # implicit import
  help_message = """Run a basic nanoBragg simulation given the imageset (detector & beam),
the dials crystal model, and the reindexed & curated strong spots.
Results can be viewed with dials.image_viewer <token>_%%%.hdf5 mask=<token>_%%%.mask
"""
  usage = ""
  from dials.util.options import OptionParser
  # Create the parser
  parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)

  # Parse the command line. quick_parse is required for MPI compatibility
  params, options = parser.parse_args(show_diff_phil=True,quick_parse=True)
  return params,options

def run(params):
    basename = "%s_%05d."%(params.output.label, params.output.index)
    M = basic_run_manager.from_files(params.trusted_mask, params.refl, params.expt)
    M.get_trusted_and_refl_mask()
    M.refl_analysis(params.cryst) # new
    M.simple_rmsd() # new
    #M.plot_pixel_histograms() # new
    M.get_lunus_repl()
    M.get_image_res_data()
    M.modify_shoeboxes() # new
    M.view["sim_mock"] = M.simulation_mockup(M.view["exp_data"]) # new
    nanobragg_sim = M.ersatz_MCMC(params.model) # initial Bragg simulation
    M.view["bragg_plus_background"] = M.reusable_rmsd(proposal=nanobragg_sim, label="ersatz_mcmc")
    M.view["renormalize_bragg_plus_background"] = M.reusable_rmsd(proposal=M.renormalize(
            proposal=nanobragg_sim,proposal_label="ersatz_mcmc",ref_label="spots_mockup"),
            label="renormalize_mcmc")
    M.view["Z_plot"] = M.Z_statistics(experiment=M.view["sim_mock"],
                                           model=M.view["renormalize_bragg_plus_background"],
                                   plot=False)
    M.write_hdf5(os.path.join(params.output.output_dir,basename+"hdf5"))
    M.resultant_mask_to_file(os.path.join(params.output.output_dir,basename+"mask"))

if __name__ == "__main__":
  params,options = parse_input()
  run(params)

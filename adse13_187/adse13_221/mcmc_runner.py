from __future__ import division
import os, pickle
from dials.array_family import flex
from scipy import constants
import copy
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt

from LS49.adse13_187.adse13_221.basic_runner import basic_run_manager
class mcmc_run_manager(basic_run_manager):

  def quick_Zscore(self, kernel_model, ref_label="spots_mockup", plot=True):
    keV_per_photon = ENERGY_CONV/1000./self.expt.beam.get_wavelength()

    # do NOT assume the model data represent the full detector image: ipanel x islow x ifast,
    # instead they are the 1D preselected active pixels from the full image, in specified order
    self.pixel_stats.analyze3(whitelist_pixels = kernel_model,
                              reference_shoebox_sums = self.refl_table[ref_label+"_shoebox_sum"],
                              slow_size = 254, panel_size = 254 * 254,
                              keV_per_photon = keV_per_photon)
    LLG = self.pixel_stats.get_LLG()
    mnz = self.pixel_stats.get_mnz()
    sgz = self.pixel_stats.get_sgz()
    proposal_center_of_mass = self.pixel_stats.get_proposal_center_of_mass()
    self.refl_table["c_temp_values"] = proposal_center_of_mass
    rmsd = self.simple_rmsd(calc_data="c_temp_values",plot=False)
    return rmsd, sgz, LLG

  def ersatz_MCMC(self, params):
    relevant_whitelist_order = flex.size_t()
    for sidx in range(len(self.refl_table["spots_offset"])):#loop through the shoeboxes
      for pidx in range(self.refl_table["spots_offset"][sidx],
                        self.refl_table["spots_offset"][sidx]+self.refl_table["spots_size"][sidx]):
        relevant_whitelist_order.append(self.spots_pixels[pidx])

    from LS49.adse13_187.adse13_221.basic_runner import MCMC_manager
    from LS49.adse13_187.adse13_221.case_run import case_job_runner
    from LS49.adse13_187.adse13_221.case_chain import case_chain_runner

    modality = "chain"

    class implement_mcmc(MCMC_manager, case_chain_runner): pass
    self.impl_MCMC = implement_mcmc()
    self.impl_MCMC.get_amplitudes(self.dials_model, self.refl_table)
    self.impl_MCMC.set_whitelist(relevant_whitelist_order)
    self.whitelist_lunus_filtered_data = []
    self.whitelist_exp_data = []
    self.whitelist_sim_mock = []
    for sidx in range(len(self.refl_table)): #loop through the shoeboxes
      for ipanel, islow, ifast in self.per_shoebox_whitelist_iterator(sidx):
        self.whitelist_lunus_filtered_data.append(self.view["lunus_filtered_data"][ipanel,islow,ifast])
        self.whitelist_exp_data.append(float(self.view["exp_data"][ipanel][islow,ifast]))
        self.whitelist_sim_mock.append(self.view["sim_mock"][ipanel][islow,ifast])
    from simtbx.pixel import pixel_stats
    self.pixel_stats = pixel_stats()
    self.pixel_stats.set_whitelist(lunus_filtered_data = self.whitelist_lunus_filtered_data,
                                   exp_data = self.whitelist_exp_data,
                                   sim_mock = self.whitelist_sim_mock)
    # proved that the above lists are defined
    self.pixel_stats.set_shoebox_iterator(shoebox_offset=self.refl_table["spots_offset"],
                                          shoebox_size=self.refl_table["spots_size"],
                                          spots_pixels=self.spots_pixels)

    C = self.impl_MCMC.chain_runner(expt=self.expt, alt_expt=self.dials_model, params=params.model,
      mask_array = self.monolithic_mask_whole_detector_as_1D_bool,
      n_cycles = params.mcmc.cycles, s_cycles = params.simplex.cycles,
      Zscore_callback=self.quick_Zscore,
      ) # do not return yet; simulated image as numpy array

    if params.mcmc.modality == "chain":
      # 1) create copy of self.dials_model, to contain updated crystal parameters from the MCMC simulation C
      bridge_model = copy.deepcopy(self.dials_model)

      # cell lengths a,c
      STM = self.impl_MCMC.parameters["cell"].get_stationary_crystal_model(self.dials_model.crystal)
          #print(STM.get_unit_cell())
          #print(self.dials_model.crystal.get_unit_cell())
      # orientation, perturbations away from dials model
      STM2 = self.impl_MCMC.parameters2["rot"].get_stationary_crystal_model(STM)
      bridge_model.crystal = STM2
          #STM2.show()
          #self.dials_model.crystal.show()

      # Ncells abc
      params.model.Nabc.value = self.impl_MCMC.parameters2["ncells"].get_stationary_model()

      # Eta abc
      mos_aniso = [self.impl_MCMC.parameters[ang].get_stationary_model() for ang in ["etaa","etab","etac"]]

      # 2) write this expt out to a file
      bridge_model.Nabc = params.model.Nabc.value
      bridge_model.mos_aniso = mos_aniso # FIXME later these properties will be part of child-class crystal model
      basename = "%s_%05d."%(params.output.label, params.output.index)
      outfile = os.path.join(params.output.output_dir,basename+"pickle")
      with open(outfile,"wb") as F: pickle.dump(bridge_model, F)

      # 3) pickle up C.parameters for offline analysis (FIXME not yet implemented)

    class ersatz(MCMC_manager, case_job_runner): pass
    self.MCMC = ersatz()
    self.MCMC.get_amplitudes(self.dials_model, self.refl_table)
    self.MCMC.set_whitelist(relevant_whitelist_order)

    if params.mcmc.modality == "job":
      return self.MCMC.job_runner(expt=self.expt, alt_expt=self.dials_model, params=params.model,
      mask_array = self.monolithic_mask_whole_detector_as_1D_bool,
      ) # returns simulated image as numpy array
    else:
      assert params.mcmc.modality == "chain"
      return self.MCMC.job_runner(expt=self.expt, alt_expt=bridge_model, params=params.model,
      mask_array = self.monolithic_mask_whole_detector_as_1D_bool,
      mos_aniso = mos_aniso) # returns simulated image as numpy array

def generate_phil_scope():
  from iotbx.phil import parse
  master_phil="""
    include scope LS49.adse13_187.adse13_221.basic_runner.phil_scope
    mcmc {
      cycles = 1000
        .type = int (value_min=10)
        .help = total number of cycles allowed for Monte Carlo
      modality = *chain job
        .type = choice
        .help = expert only. Chain means finalize the MCMC results, publish expt, output pickle, write
        .help = hdf5, or job means do the analysis for the input crystal model, not the MCMC results.
    }
    simplex {
      cycles = 0
        .type = int (value_min=0)
        .help = maximum simplex cycles to iterate prior to mcmc
    }
    model {
      cell {
        covariance = ""
          .type = path
          .help = pickle file of the covariance model as described in the tdata documentation for cctbx.xfel.merge
        hyperparameter = 0.5
          .type = float(value_min=0.01, value_max=5)
          .help = delta-value proposal for the next markov chain iteration, expressed in covariance matrix
          .help = Mahalanobis units (number of standard deviations to scale the jump).
      }
      mosaic_spread {
        hyperparameter = 0.2
          .type = float(value_min=0.001, value_max=0.9)
          .help = The allowable range of proposal values as a multiplier of current value.
          .help = If the current value is 1.0, a hyperparameter of 0.2 would give an allowable range of 0.8-1.2.
      }
      rot {
        refine = False
          .type = bool
        hyperparameter = 0.04
          .type = float(value_min=0.001, value_max=0.9)
          .help = The allowable range of proposal values as a multiplier of current value.
        sigmas = (0.03,0.03,0.01)
        .type = floats(size=3, value_min=0)
        .help = Starting allowable range in degrees for rotx,roty,rotz
        value = (0.00,0.00,0.00)
        .type = floats(size=3)
        .help = Starting allowable range in degrees for rotx,roty,rotz
      }
    }
  """
  return parse(master_phil, process_includes=True)
phil_scope = generate_phil_scope()

def parse_input():
  # The script usage
  import libtbx.load_env # implicit import
  help_message = """Upgrade the basic nanoBragg simulation to an Markov Chain Monte Carlo
run where Stage 1 parameters are fit.
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
    M = mcmc_run_manager.from_files(params.trusted_mask, params.refl, params.expt)
    M.remove_overlaps()
    M.get_trusted_and_refl_mask()
    M.refl_analysis(params.cryst) # new, sets M.dials_model from params.cryst
    M.simple_rmsd() # new
    #M.plot_pixel_histograms() # new
    M.get_lunus_repl()
    M.get_image_res_data()
    M.modify_shoeboxes() # new
    M.view["sim_mock"] = M.simulation_mockup(M.view["exp_data"]) # new
    nanobragg_sim = M.ersatz_MCMC(params = params) # final Bragg simulation after MCMC run
    M.view["bragg_plus_background"] = M.reusable_rmsd(proposal=nanobragg_sim, label="ersatz_mcmc")
    M.view["renormalize_bragg_plus_background"] = M.reusable_rmsd(proposal=M.renormalize(
            proposal=nanobragg_sim,proposal_label="ersatz_mcmc",ref_label="spots_mockup"),
            label="renormalize_mcmc",plot=True)
    M.view["Z_plot"] = M.Z_statistics(experiment=M.view["sim_mock"],
                                           model=M.view["renormalize_bragg_plus_background"],
                                   plot=True)
    M.write_hdf5(os.path.join(params.output.output_dir,basename+"hdf5"))
    M.resultant_mask_to_file(os.path.join(params.output.output_dir,basename+"mask"))

if __name__ == "__main__":
  params,options = parse_input()
  run(params)

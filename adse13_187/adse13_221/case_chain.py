from __future__ import division
from LS49.adse13_187.cyto_batch import multipanel_sim
from time import time
import random

class case_chain_runner:
  def chain_runner(self,expt,mask_array=None,n_cycles = 100,
      Zscore_callback=None, rmsd_callback=None):

    # Fixed hyperparameters
    mosaic_spread_samples = 250
    beamsize_mm = 0.000886226925452758 # sqrt beam focal area
    spot_scale = 500.
    oversample = 1  # factor 1,2, or 3 probably enough
    verbose = 0  # leave as 0, unless debug
    shapetype = "gauss_argchk"

    detector = expt.detector
    flat = True  # enforce that the camera has 0 thickness
    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)
        assert detector[0].get_thickness() == 0

    beam = expt.beam
    spec = expt.imageset.get_spectrum(0)
    energies_raw = spec.get_energies_eV().as_numpy_array()
    weights_raw = spec.get_weights().as_numpy_array()
    from LS49.adse13_187.adse13_221.explore_spectrum import method3
    energies, weights, _ = method3(energies_raw, weights_raw,); weights = 5000000.*weights
    energies = list(energies); weights = list(weights)

    device_Id = 0 # XXX revisit for multiprocess service
    assert self.gpu_channels_singleton is not None
    device_Id = self.gpu_channels_singleton.get_deviceID()
    print (device_Id,"device",shapetype)
    Famp_is_uninitialized = ( self.gpu_channels_singleton.get_nchannels() == 0 )
    if Famp_is_uninitialized:
      F_P1 = self.amplitudes
      for x in range(1):  # in this scenario, amplitudes are independent of lambda
          self.gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
    assert self.gpu_channels_singleton.get_nchannels() == 1

    # Variable parameters
    mosaic_spread = 0.01  # degrees
    Ncells_abc = 130, 30, 10  # medians from best stage1
    alt_exper = ExperimentListFactory.from_json_file(
      '/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt',
                                              check_format=False)[0]
    alt_crystal = alt_exper.crystal
    from LS49.adse13_187.adse13_221.parameters import variable_cell, variable_mosaicity, covariant_cell
    self.parameters = {}
    #self.parameters["cell"] = variable_cell(alt_crystal)
    self.parameters["cell"] = covariant_cell.from_covariance(alt_crystal)
    #self.parameters["eta"] = variable_mosaicity(mosaic_spread)
    self.parameters["etaa"] = variable_mosaicity(mosaic_spread)
    self.parameters["etab"] = variable_mosaicity(mosaic_spread)
    self.parameters["etac"] = variable_mosaicity(mosaic_spread)
    self.rmsd_chain= flex.double(); self.sigz_chain=flex.double(); self.llg_chain=flex.double();
    self.cycle_list = [key for key in self.parameters]
    self.accept = flex.int()

    for macro_iteration in range(n_cycles):
      BEG=time()
      turn = self.cycle_list[macro_iteration%len(self.cycle_list)]
      if turn=="cell":
        alt_crystal = self.parameters["cell"].get_current_crystal_model()

      whitelist_only, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=energies, fluxes=weights,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=self.parameters["etaa"].proposal,
        mos_aniso=(self.parameters["etaa"].proposal,0,0,0,self.parameters["etab"].proposal,0,0,0,self.parameters["etac"].proposal),
        beamsize_mm=beamsize_mm,
        profile=shapetype,
        show_params=False,
        time_panels=False, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=False,
        mask_file=mask_array, skip_numpy=True,
        relevant_whitelist_order=self.relevant_whitelist_order
      )
      Rmsd,sigZ,LLG = Zscore_callback(kernel_model=whitelist_only, plot=False)
      if macro_iteration==0:
        for key in self.parameters:
          self.parameters[key].accept()
        self.accept.append(1)
        self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)
      else:
        print ("Old NLL ",self.llg_chain[-1], "NEW LLG",LLG, "diff",self.llg_chain[-1] - LLG)
        this_cycle_key = self.cycle_list[(macro_iteration)%len(self.cycle_list)]
        acceptance_prob = min (
          1.,
          math.exp( (self.llg_chain[-1] - LLG)/55000. ) # normalize by number of pixels # XXX FIXME denominator
          * self.parameters[this_cycle_key].transition_probability_ratio # q(X|Y)/q(Y|X), Y=proposal, X=last value
        )

        if random.random() < acceptance_prob:
          for key in self.parameters:
            if key == turn:  self.parameters[key].accept()
            else: self.parameters[key].reject()
          self.accept.append(1)
          self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)
        else:
          for key in self.parameters: self.parameters[key].reject()
          self.accept.append(0)
          self.rmsd_chain.append(self.rmsd_chain[-1]);
          self.sigz_chain.append(self.sigz_chain[-1]);
          self.llg_chain.append(self.llg_chain[-1])

      for key in self.parameters:
        if key == self.cycle_list[(macro_iteration+1)%len(self.cycle_list)]:
          self.parameters[key].generate_next_proposal()
      self.plot_all(macro_iteration+1,of=n_cycles)
      TIME_EXA = time()-BEG
      print("\t\tExascale: time for Bragg sim: %.4fs; total: %.4fs" % (TIME_BRAGG, TIME_EXA))
    print ("Mean RMSD %.2f"%(flex.mean(self.rmsd_chain[len(self.rmsd_chain)//2:])))
    print ("Mean sigz %.2f"%(flex.mean(self.sigz_chain[len(self.sigz_chain)//2:])))
    print ("Mean -LLG %.2f"%(flex.mean(self.llg_chain[len(self.llg_chain)//2:])))
    exit("XXX still need to implement return values from MCMC")
    return JF16M_numpy_array.as_numpy_array()


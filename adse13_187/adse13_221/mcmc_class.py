from __future__ import division
import os
from time import time
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from LS49.adse13_187.cyto_batch import multipanel_sim
from matplotlib import pyplot as plt
import random, math
from scipy import constants
ENERGY_CONV = 1e10*constants.c*constants.h / constants.electron_volt

class MCMC_manager:
  def get_amplitudes(self, dials_model, refl_table, test_without_mpi=True):
    from LS49.adse13_187.cyto_batch import parse_input
    self.params,options = parse_input()

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

  def chain_runner(self,expt,mask_array=None,n_cycles = 100,
      Zscore_callback=None, rmsd_callback=None):

    # Fixed hyperparameters
    mosaic_spread_samples = 250
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    verbose = 0  # leave as 0, unles debug
    shapetype = "gauss_argchk"

    detector = expt.detector
    flat = True  # enfore that the camera has 0 thickness
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

    device_Id = 0 # XXX probably has to be revisited for multiprocess service
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
    self.parameters["eta"] = variable_mosaicity(mosaic_spread)
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
        mos_dom=mosaic_spread_samples, mos_spread=self.parameters["eta"].proposal,
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
        self.parameters["cell"].accept()
        self.parameters["eta"].accept()
        self.accept.append(1)
        self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)
      else:
        print ("Old NLL ",self.llg_chain[-1], "NEW LLG",LLG, "diff",self.llg_chain[-1] - LLG)
        acceptance_prob = min (1., math.exp( (self.llg_chain[-1] - LLG)/55000. )) # normalize by number of pixels

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
    exit("XXX still need to implement return values from MCMC")
    return JF16M_numpy_array.as_numpy_array()

  def plot_all(self,icmp,of):
    N_param = len(self.parameters)
    if icmp==1:
      self.line = list(range(N_param))
      #self.running = flex.double(); self.running.append(self.accept[0])
      plt.ion()
      self.fig,self.axes = fig,axes = plt.subplots(7,1,sharex=True,figsize=(7,10))
      self.line0, = axes[0].plot(range(icmp), self.parameters["cell"].a_chain)
      self.line1, = axes[1].plot(range(icmp), self.parameters["cell"].c_chain)
      self.line1eta, = axes[2].plot(range(icmp), self.parameters["cell"].c_chain)
      self.line2, = axes[3].plot(range(icmp), self.rmsd_chain)
      self.line3, = axes[4].plot(range(icmp), self.sigz_chain)
      self.line4, = axes[5].plot(range(icmp), self.llg_chain)
      #self.line5, = axes[6].plot(range(icmp), self.running, "r-")
      self.line5a, = axes[6].plot(range(icmp), self.accept, "k,")
      for npm in range(N_param):
        self.line[npm], = axes[6].plot(range(icmp), N_param*self.parameters[self.cycle_list[npm]].running)
      print(self.line)
      axes[0].set_ylabel("a")
      axes[1].set_ylabel("c")
      axes[2].set_ylabel("eta")
      axes[3].set_ylabel("RMSD (px)")
      axes[4].set_ylabel("Z-score sigma")
      axes[5].set_ylabel("LLG")
      axes[6].set_ylabel("accept")
      axes[0].set_ylim(78.4,79.1)
      axes[1].set_ylim(264.5,266)
      axes[2].set_ylim(0,0.1)
      axes[3].set_ylim(0.45,1.00)
      axes[4].set_ylim(2.0,6.0)
      axes[5].set_ylim(-17052536,-16052536)
      axes[6].set_ylim(-0.1,1.1)
      plt.xlim(0,of+1)
      plt.show()
    elif icmp%10==0:
#start here
#1) use multivariate gaussian
#2) add rotations
#3) take some easy way of writing out the results & generating a report.
#4) look over the callback for any quick way to double the speed, such as pre-generating lists
      #self.running.append(flex.sum(self.accept)/len(self.accept))
      self.line0.set_xdata(range(icmp))
      self.line0.set_ydata(self.parameters["cell"].a_chain)
      self.line1.set_xdata(range(icmp))
      self.line1.set_ydata(self.parameters["cell"].c_chain)
      self.line1eta.set_xdata(range(icmp))
      self.line1eta.set_ydata(self.parameters["eta"].chain)
      self.line2.set_xdata(range(icmp))
      self.line2.set_ydata(self.rmsd_chain)
      self.line3.set_xdata(range(icmp))
      self.line3.set_ydata(self.sigz_chain)
      self.line4.set_xdata(range(icmp))
      self.line4.set_ydata(self.llg_chain)
      #self.line5.set_xdata(range(icmp))
      #self.line5.set_ydata(self.running)
      self.line5a.set_xdata(range(icmp))
      self.line5a.set_ydata(self.accept)
      for npm in range(N_param):
        self.line[npm].set_xdata(range(icmp))
        self.line[npm].set_ydata(N_param*self.parameters[self.cycle_list[npm]].running)
      self.fig.canvas.draw()
      self.fig.canvas.flush_events()
      if icmp==of-1: input()

  def job_runner(self,expt,mask_array=None,i_exp=0,spectra={}):

    mosaic_spread = 0.00  # degrees
    mosaic_spread_samples = 250 # XXX Fixme make this a parameter
    Ncells_abc = 130, 30, 10  # medians from best stage1
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    #<><><><><><><><>
    exper = expt

    crystal = exper.crystal
    detector = exper.detector

    alt_exper = ExperimentListFactory.from_json_file(
      '/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt',
                                              check_format=False)[0]
    AC = alt_crystal = alt_exper.crystal

    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)

    beam = exper.beam

    spec = exper.imageset.get_spectrum(0)
    energies_raw, weights_raw = spec.get_energies_eV().as_numpy_array(), \
                                spec.get_weights().as_numpy_array()
    from LS49.adse13_187.adse13_221.explore_spectrum import method3
    energies, weights, _ = method3(energies_raw, weights_raw,); weights = 5000000.*weights

    if flat:
        assert detector[0].get_thickness() == 0

    device_Id = 0
    if self.gpu_channels_singleton is not None:
      device_Id = self.gpu_channels_singleton.get_deviceID()

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = ENERGY_CONV / mn_energy
    print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("\tBreakdown:")
    for shapetype in ["gauss_argchk"]:
      BEG=time()
      print (self.gpu_channels_singleton.get_deviceID(),"device",shapetype)
      Famp_is_uninitialized = ( self.gpu_channels_singleton.get_nchannels() == 0 )
      if Famp_is_uninitialized:
        F_P1 = self.amplitudes
        for x in range(1):  # in this scenario, amplitudes are independent of lambda
          self.gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
      assert self.gpu_channels_singleton.get_nchannels() == 1

      from LS49.adse13_187.cyto_batch import multipanel_sim
      JF16M_numpy_array, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=list(energies), fluxes=list(weights),
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
        beamsize_mm=beamsize_mm,
        profile=shapetype,
        show_params=False,
        time_panels=False, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=False,
        mask_file=(mask_array or self.params.mask_file))
      TIME_EXA = time()-BEG

      print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    return JF16M_numpy_array

# use case specializations are given below:
class case_DS1 (MCMC_manager):
  # calculate an rmsd on these pixels (omit background) DONE
  # try recruiting Derek's March 15 expt model to see if rmsd improves DONE
  # renormalize the proposal for each Bragg spot to be on scale with experiment DONE
  # calculate a Z-score, display the image, also calculate a per-spot and all-spot probability DONE
  # remove background calculation DONE

  def job_runner(self,i_exp=0,spectra={}):
    from simtbx.nanoBragg import utils
    print("Experiment %d" % i_exp, flush=True)

    from LS49.adse13_187.case_data import retrieve_from_repo
    experiment_file = retrieve_from_repo(i_exp)
    cuda = True  # False  # whether to use cuda
    mosaic_spread = 0.00  # degrees
    mosaic_spread_samples = 250 # XXX Fixme make this a parameter
    Ncells_abc = 130, 30, 10  # medians from best stage1
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    include_background = False
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    #<><><><><><><><>
    os.environ["NXMX_LOCAL_DATA"]="/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5"
    El = ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=True)
    exper = El[0]

    crystal = exper.crystal
    detector = exper.detector

    alt_exper = ExperimentListFactory.from_json_file(
      '/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt',
                                              check_format=False)[0]
    AC = alt_crystal = alt_exper.crystal

    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)

    beam = exper.beam

    # XXX new code
    spec = exper.imageset.get_spectrum(0)
    energies_raw, weights_raw = spec.get_energies_eV().as_numpy_array(), \
                                spec.get_weights().as_numpy_array()
    energies, weights = utils.downsample_spectrum(energies_raw, weights_raw, method=1, total_flux=total_flux,
                                                  ev_width=ev_res)

    if flat:
        assert detector[0].get_thickness() == 0

    device_Id = 0
    if self.gpu_channels_singleton is not None:
      device_Id = self.gpu_channels_singleton.get_deviceID()

    show_params = False

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = utils.ENERGY_CONV / mn_energy
    print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("\tBreakdown:")
    for shapetype in ["gauss_argchk"]:
      BEG=time()
      print (self.gpu_channels_singleton.get_deviceID(),"device",shapetype)
      Famp_is_uninitialized = ( self.gpu_channels_singleton.get_nchannels() == 0 )
      if Famp_is_uninitialized:
        F_P1 = self.amplitudes
        for x in range(1):  # in this scenario, amplitudes are independent of lambda
          self.gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
      assert self.gpu_channels_singleton.get_nchannels() == 1

      from LS49.adse13_187.cyto_batch import multipanel_sim
      JF16M_numpy_array, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=list(energies), fluxes=list(weights),
        background_wavelengths=[mn_wave], background_wavelength_weights=[1],
        background_total_flux=total_flux,background_sample_thick_mm=0.5,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
        beamsize_mm=beamsize_mm,
        profile=shapetype,
        show_params=show_params,
        time_panels=False, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=include_background,
        mask_file=self.params.mask_file)
      TIME_EXA = time()-BEG

      print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    return JF16M_numpy_array

class case_228 (MCMC_manager):
  def job_runner(self,i_exp=0,spectra={}):
    from simtbx.nanoBragg import utils
    print("Experiment %d" % i_exp, flush=True)

    from LS49.adse13_187.case_data import retrieve_from_repo
    experiment_file = retrieve_from_repo(i_exp)
    cuda = True  # False  # whether to use cuda
    mosaic_spread = 0.07  # degrees
    mosaic_spread_samples = 500 # XXX Fixme make this a parameter
    Ncells_abc = 30, 30, 10  # medians from best stage1
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    include_background = False
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    #<><><><><><><><>
    os.environ["NXMX_LOCAL_DATA"]="/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5"
    El = ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=True)
    exper = El[0]

    crystal = exper.crystal
    detector = exper.detector
    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)

    beam = exper.beam

    # XXX new code
    spec = exper.imageset.get_spectrum(0)
    energies_raw, weights_raw = spec.get_energies_eV().as_numpy_array(), \
                                spec.get_weights().as_numpy_array()
    energies, weights = utils.downsample_spectrum(energies_raw, weights_raw, method=1, total_flux=total_flux,
                                                  ev_width=ev_res)

    if flat:
        assert detector[0].get_thickness() == 0

    device_Id = 0
    if self.gpu_channels_singleton is not None:
      device_Id = self.gpu_channels_singleton.get_deviceID()

    show_params = False

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = utils.ENERGY_CONV / mn_energy
    print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("\tBreakdown:")
    for shapetype in ["gauss_argchk"]:
      BEG=time()
      print (self.gpu_channels_singleton.get_deviceID(),"device",shapetype)
      Famp_is_uninitialized = ( self.gpu_channels_singleton.get_nchannels() == 0 )
      if Famp_is_uninitialized:
        from iotbx.reflection_file_reader import any_reflection_file
        from LS49 import ls49_big_data
        merge_file = os.path.join(ls49_big_data,"adse13_228","cyto_init_merge.mtz")
        self.merged_amplitudes = any_reflection_file(merge_file).as_miller_arrays()[0].as_amplitude_array()

        F1 = self.merged_amplitudes.expand_to_p1()
        F2 = self.amplitudes.expand_to_p1() # takes care of both transform to asu & expand

        if False: # make sure that mtz file (F1) and strong spots (self.amplitudes) are roughly correlated
          from matplotlib import pyplot as plt
          from cctbx import miller
          matches = miller.match_indices( F1.indices(), self.amplitudes.indices() )
          sel0 = flex.size_t([p[0] for p in matches.pairs()])
          sel1 = flex.size_t([p[1] for p in matches.pairs()])
          data0 = F1.data().select(sel0)
          data1 = self.amplitudes.data().select(sel1)
          plt.plot(data0, data1, 'r.')
          plt.show() # yes, the two are very roughly correlated
          # end of test

        #F_P1 = F1 # legacy, use a merged mtz file
        #F_P1 = F2 # this one way absolutely wrong! way too many predictions, beyond the strong spots
        F_P1 = F1
        for x in range(1):  # in this scenario, amplitudes are independent of lambda
          self.gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
      assert self.gpu_channels_singleton.get_nchannels() == 1

      from LS49.adse13_187.cyto_batch import multipanel_sim
      JF16M_numpy_array, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=list(energies), fluxes=list(weights),
        background_wavelengths=[mn_wave], background_wavelength_weights=[1],
        background_total_flux=total_flux,background_sample_thick_mm=0.5,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
        beamsize_mm=beamsize_mm,
        profile=shapetype,
        show_params=show_params,
        time_panels=False, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=include_background,
        mask_file=self.params.mask_file)
      TIME_EXA = time()-BEG

      print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    return JF16M_numpy_array

def proof_of_principle_compare_three_spectra(energies_raw, weights_raw):
    from simtbx.nanoBragg import utils
    energies1, weights1 = utils.downsample_spectrum(energies_raw, weights_raw, method=1, total_flux=1e12,
                                                  ev_width=1.5)
    energies2, weights2 = utils.downsample_spectrum(energies_raw, weights_raw, method=2, total_flux=1e12,
                                                  ev_width=1.5)
    from LS49.adse13_187.adse13_221.explore_spectrum import method3
    energies, weights, _ = method3(energies_raw, weights_raw,); weights = 5000000.*weights

    from matplotlib import pyplot as plt
    plt.plot(energies1,weights1,'r-')
    plt.plot(energies2,weights2,'b-')
    plt.plot(energies, weights, 'g-')
    plt.show()
    return energies, weights


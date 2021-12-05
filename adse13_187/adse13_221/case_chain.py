from __future__ import division
from LS49.adse13_187.cyto_batch import multipanel_sim
from time import time
import random, math
from scitbx.array_family import flex
from matplotlib import pyplot as plt
import copy
from scitbx.matrix import sqr, col
cube_diag = math.sqrt(1./3) # 0.57735
unit_vectors = [(1,0,0), (0,1,0), (0,0,1), (cube_diag, cube_diag, cube_diag)]

class case_chain_runner:
  def chain_runner(self,expt,alt_expt,params,mask_array=None,n_cycles = 100, s_cycles=0,
      Zscore_callback=None, rmsd_callback=None):

    self.model_plot_enable = params.plot
    if mask_array is not None:
      assert type(mask_array) is flex.bool # type check intending to convert active-pixel-bools to whitelist-ints
      active_pixels = flex.int()
      for i, x in enumerate(mask_array):
        if x: active_pixels.append(i)
    mask_array = active_pixels

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

    alt_crystal = copy.deepcopy(alt_expt.crystal) # avoid perturbing the original dials cell

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
          self.gpu_channels_singleton.structure_factors_to_GPU_direct(
          x, F_P1.indices(), F_P1.data())
    assert self.gpu_channels_singleton.get_nchannels() == 1

    # Variable parameters
    mosaic_spread = params.mosaic_spread.value
    Ncells_abc = params.Nabc.value

    from LS49.adse13_187.adse13_221.parameters import variable_mosaicity
    from LS49.adse13_187.adse13_221.parameters import covariant_cell, covariant_rot, covariant_ncells
    self.parameters = {}
    self.parameters["cell"] = covariant_cell.from_covariance(alt_crystal, params.cell)
    self.parameters["etaa"] = variable_mosaicity(mosaic_spread, label="η a", params = params.mosaic_spread)
    self.parameters["etab"] = variable_mosaicity(mosaic_spread, label="η b", params = params.mosaic_spread)
    self.parameters["etac"] = variable_mosaicity(mosaic_spread, label="η c", params = params.mosaic_spread)
    self.parameters2 = {}
    if params.rot.refine:
      self.parameters2["rot"] = covariant_rot(alt_crystal, params.rot)
    if params.Nabc.refine:
      self.parameters2["ncells"] = covariant_ncells(params.Nabc)
    self.ref_params={};self.ref_params.update(self.parameters); self.ref_params.update(self.parameters2)
    # XXX TO DO list (Nick/Dan discuss)
    # 1) change the variable mosaicity model to use updated aniso Derek model (Nick)
    # 2) add rotx/roty/rotz.  Covariant excursion values from Sauter 2014 paper: rotz 0.02° rotx 0.03° roty 0.03°
    # 3) refine ncells a b c

    self.rmsd_chain= flex.double(); self.sigz_chain=flex.double(); self.llg_chain=flex.double();
    self.cycle_list = [key for key in self.ref_params]
    self.accept = flex.int()

    self.beginning_iteration = 0
    if s_cycles > 0:
      from LS49.adse13_187.adse13_221.simplex_method import simplex_detail
      # initialize prior to simplex
      whitelist_only, TIME_BG, TIME_BRAGG, self.exascale_mos_blocks = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=energies, fluxes=weights,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=self.parameters["etaa"].proposal,
        mos_aniso=(self.parameters["etaa"].proposal,
                   self.parameters["etab"].proposal,
                   self.parameters["etac"].proposal),
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
      self.accept.append(1)
      self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)

      PP = dict(detector=detector, beam=beam,
        energies=energies, weights=weights,
        oversample=oversample,
        mosaic_spread_samples=mosaic_spread_samples,
        beamsize_mm=beamsize_mm,
        shapetype=shapetype,
        verbose=verbose,
        spot_scale=spot_scale,
        mask_array=mask_array, Z=Zscore_callback)

      MIN = simplex_detail(alt_crystal, Ncells_abc, host_runner = self,
             PP = PP, n_cycles=n_cycles, s_cycles=s_cycles)
      self.beginning_iteration = MIN.iteration+1

    for macro_iteration in range(self.beginning_iteration, n_cycles):
      BEG=time()
      turn = self.cycle_list[macro_iteration%len(self.cycle_list)]
      if turn=="cell":
        alt_crystal = self.parameters["cell"].get_current_crystal_model(alt_crystal)
      elif turn=="rot":
        alt_crystal = self.parameters2["rot"].get_current_crystal_model(alt_crystal)
      elif turn=="ncells":
        Ncells_abc = self.parameters2["ncells"].get_current_model()

      whitelist_only, TIME_BG, TIME_BRAGG, self.exascale_mos_blocks = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=energies, fluxes=weights,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=self.parameters["etaa"].proposal,
        mos_aniso=(self.parameters["etaa"].proposal,
                   self.parameters["etab"].proposal,
                   self.parameters["etac"].proposal),
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
      if macro_iteration==self.beginning_iteration:
        for key in self.ref_params:
          self.ref_params[key].accept()
        self.accept.append(1)
        self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)
      else:
        print ("Old NLL ",self.llg_chain[-1], "NEW LLG",LLG, "diff",self.llg_chain[-1] - LLG)
        this_cycle_key = self.cycle_list[(macro_iteration)%len(self.cycle_list)]
        acceptance_prob = min (
          1.,
          math.exp( (self.llg_chain[-1] - LLG)/len(whitelist_only) ) # normalize by no. of pixels
          * self.ref_params[this_cycle_key].transition_probability_ratio # q(X|Y)/q(Y|X), Y=proposal, X=last value
        )
        if random.random() < acceptance_prob:
          for key in self.ref_params:
            if key == turn:  self.ref_params[key].accept()
            else: self.ref_params[key].reject()
          self.accept.append(1)
          self.rmsd_chain.append(Rmsd); self.sigz_chain.append(sigZ); self.llg_chain.append(LLG)
        else:
          for key in self.ref_params: self.ref_params[key].reject()
          self.accept.append(0)
          self.rmsd_chain.append(self.rmsd_chain[-1]);
          self.sigz_chain.append(self.sigz_chain[-1]);
          self.llg_chain.append(self.llg_chain[-1])
      P = Profiler("%40s"%"key maintenance")

      for key in self.ref_params:
        if key == self.cycle_list[(macro_iteration+1)%len(self.cycle_list)]:
          self.ref_params[key].generate_next_proposal()
      P = Profiler("%40s"%"plot all")
      self.plot_all(macro_iteration+1,of=n_cycles)
      del P
      TIME_EXA = time()-BEG
      print("\t\tExascale: time for Bragg sim: %.4fs; total: %.4fs\n" % (TIME_BRAGG, TIME_EXA))
    print ("MCMC <RMSD> %.2f"%(flex.mean(self.rmsd_chain[len(self.rmsd_chain)//2:])))
    print ("MCMC <sigz> %.2f"%(flex.mean(self.sigz_chain[len(self.sigz_chain)//2:])))
    print ("MCMC <-LLG> %.2f"%(flex.mean(self.llg_chain[len(self.llg_chain)//2:])))
    for key in self.ref_params: self.ref_params[key].show()
    if self.model_plot_enable:
      plt.close(self.fig)
      plt.close(self.fig2)
      plt.close(self.fig3)
      plt.ioff()
    print(flush=True)

  def plot_all(self,icmp,of):
    if not self.model_plot_enable: return
    N_param = len(self.parameters)
    if icmp==1:
      self.lines = {}
      self.line = list(range(N_param))
      ax_ct_1 = 4
      for key in self.parameters: ax_ct_1 += self.parameters[key].display_n
      plt.ion()
      self.fig,self.axes = fig,axes = plt.subplots(ax_ct_1,1,sharex=True,figsize=(7,10))
      self.line2, = axes[0].plot(range(icmp), self.rmsd_chain)
      self.line3, = axes[1].plot(range(icmp), self.sigz_chain)
      self.line4, = axes[2].plot(range(icmp), self.llg_chain)
      self.line5a, = axes[3].plot(range(icmp), self.accept, "k,")
      for npm in range(N_param):
        self.line[npm], = axes[3].plot(range(self.beginning_iteration, icmp), N_param*self.parameters[self.cycle_list[npm]].running)
      print("LINE",self.line)
      axes[0].set_ylabel("RMSD (px)")
      axes[1].set_ylabel("Z-score\nsigma")
      axes[2].set_ylabel("LLG")
      axes[3].set_ylabel("accept")
      axes[0].set_ylim(0.45,1.00)
      axes[1].set_ylim(1.0,6.0)
      axes[2].set_ylim(-17052536,-16052536) #clearly needs to be customized for each run
      axes[3].set_ylim(-0.1,1.1)

      axes_idx = 3
      for key in self.parameters:
        for label in self.parameters[key].display_labels:
          axes_idx += 1
          self.lines[key+label] = (axes[axes_idx].plot(range(icmp), self.parameters[key].chain[label]))[0]
          axes[axes_idx].set_ylabel(self.parameters[key].formatt%label)
          axes[axes_idx].set_ylim(self.parameters[key].display_ranges[label])
      plt.xlim(0,of+1)
      plt.tight_layout()
      if len(self.parameters2)>0:
        self.lines2 = {}
        ax_ct_2 = 0
        for key in self.parameters2: ax_ct_2 += self.parameters2[key].display_n

        self.fig2,self.axes2 = plt.subplots(ax_ct_2,1,sharex=True,figsize=(7,10))

        axes_idx = -1
        for key in self.parameters2:
          for label in self.parameters2[key].display_labels:
            axes_idx += 1
            self.lines2[key+label] = (self.axes2[axes_idx].plot(range(icmp), self.parameters2[key].chain[label]))[0]
            self.axes2[axes_idx].set_ylabel(self.parameters2[key].formatt%label)
            self.axes2[axes_idx].set_ylim(self.parameters2[key].display_ranges[label])
      plt.xlim(0,of+1)
      if True: #mosaic plot
        self.lines3 = {}
        self.fig3,self.axes3 = plt.subplots(1,4,sharey=True,figsize=(15,7))
        self.axes3[0].set_ylabel("unit projection")
        for icol, RLP in enumerate(unit_vectors):
          RLP = col(RLP)
          axis = self.axes3[icol]
          unit = RLP.normalize()
          seed = col(unit_vectors[(icol+2)%(len(unit_vectors))])
          perm2 = unit.cross(seed)
          perm3 = unit.cross(perm2)
          a2 = flex.double(); a3 = flex.double()
          for u in self.exascale_mos_blocks:
            U = sqr(u)
            newvec = U * unit
            a2.append(newvec.dot(perm2)); a3.append(newvec.dot(perm3))
          self.lines3[icol] = axis.plot (a2,a3,'r,')[0]
          axis.set_aspect("equal")
          axis_str = "(%5.3f,%5.3f,%5.3f)"%(RLP.elems) if icol==3 else "(%.0f,%.0f,%.0f)"%(RLP.elems)
          axis.set_title("Transformation of unit vector\n%s"%(axis_str))
          axis.set_xlabel("unit projection")
          axis.set_xlim(-0.0005,0.0005)
          axis.set_ylim(-0.0005,0.0005)
      plt.show()
    elif icmp%50==0:
      for key in self.parameters:
        for label in self.parameters[key].display_labels:
          self.lines[key+label].set_xdata(range(icmp))
          self.lines[key+label].set_ydata(self.parameters[key].chain[label])

      self.line2.set_xdata(range(icmp))
      self.line2.set_ydata(self.rmsd_chain)
      self.line3.set_xdata(range(icmp))
      self.line3.set_ydata(self.sigz_chain)
      self.line4.set_xdata(range(icmp))
      self.line4.set_ydata(self.llg_chain)
      self.line5a.set_xdata(range(icmp))
      self.line5a.set_ydata(self.accept)
      for npm in range(N_param):
        self.line[npm].set_xdata(range(icmp))
        self.line[npm].set_ydata(N_param*self.parameters[self.cycle_list[npm]].running)

      if len(self.parameters2)>0:
        for key in self.parameters2:
          for label in self.parameters2[key].display_labels:
            self.lines2[key+label].set_xdata(range(icmp))
            self.lines2[key+label].set_ydata(self.parameters2[key].chain[label])
      if True: #mosaic plot
        for icol, RLP in enumerate(unit_vectors):
          RLP = col(RLP)
          axis = self.axes3[icol]
          unit = RLP.normalize()
          seed = col(unit_vectors[(icol+2)%(len(unit_vectors))])
          perm2 = unit.cross(seed)
          perm3 = unit.cross(perm2)
          a2 = flex.double(); a3 = flex.double()
          for u in self.exascale_mos_blocks:
            U = sqr(u)
            newvec = U * unit
            a2.append(newvec.dot(perm2)); a3.append(newvec.dot(perm3))
          self.lines3[icol].set_xdata(a2)
          self.lines3[icol].set_ydata(a3)

      self.fig.canvas.draw()
      self.fig.canvas.flush_events()
      if icmp==of-1: input()


from __future__ import division
from LS49.adse13_187.cyto_batch import multipanel_sim
from time import time
from dxtbx.model.experiment_list import ExperimentListFactory
import os

class case_DS1:
  def job_runner(self,i_exp=0,spectra={}):
    from simtbx.nanoBragg import utils

    from LS49.adse13_187.case_data import retrieve_from_repo
    experiment_file = retrieve_from_repo(i_exp)

    # Fixed hyperparameters
    mosaic_spread_samples = 250
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758 # sqrt beam focal area
    spot_scale = 500.
    oversample = 1  # factor 1,2, or 3 probably enough
    verbose = 0  # leave as 0, unless debug
    shapetype = "gauss_argchk"

    #<><><><><><><><>
    os.environ["NXMX_LOCAL_DATA"]="/global/cfs/cdirs/m3562/der/master_files/run_000795.JF07T32V01_master.h5"
    expt = ExperimentListFactory.from_json_file(
      experiment_file, check_format=True)[0]

    crystal = expt.crystal
    detector = expt.detector
    flat = True  # enforce that the camera has 0 thickness
    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)
        assert detector[0].get_thickness() == 0

    alt_exper = ExperimentListFactory.from_json_file(
      '/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt',
      check_format=False)[0]
    AC = alt_crystal = alt_exper.crystal

    beam = expt.beam
    spec = expt.imageset.get_spectrum(0)
    energies_raw = spec.get_energies_eV().as_numpy_array()
    weights_raw = spec.get_weights().as_numpy_array()
    energies, weights = utils.downsample_spectrum(
       energies_raw, weights_raw, method=1, total_flux=total_flux, ev_width=ev_res)

    device_Id = 0
    if self.gpu_channels_singleton is not None:
      device_Id = self.gpu_channels_singleton.get_deviceID()

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

      # Variable parameters
      mosaic_spread = 0.00 # degrees
      Ncells_abc = 130, 30, 10  # medians from best stage1

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
        show_params=False,
        time_panels=False, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=False,
        mask_file=mask_array)
      TIME_EXA = time()-BEG

      print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
    return JF16M_numpy_array

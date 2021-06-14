from __future__ import division
from LS49.adse13_187.cyto_batch import multipanel_sim
from time import time
from scitbx.array_family import flex

class case_job_runner:
  def job_runner(self,expt,alt_expt,params,mask_array=None,i_exp=0,spectra={}):

    # Fixed hyperparameters
    mosaic_spread_samples = 250
    beamsize_mm = 0.000886226925452758 # sqrt beam focal area
    spot_scale = 500.
    oversample = 1  # factor 1,2, or 3 probably enough
    verbose = 0  # leave as 0, unless debug
    shapetype = "gauss_argchk"

    if mask_array is not None:
      active_pixels = flex.int()
      for i, x in enumerate(mask_array):
        if x: active_pixels.append(i)
      mask_array = active_pixels

    detector = expt.detector
    flat = True  # enforce that the camera has 0 thickness
    if flat:
        from dxtbx_model_ext import SimplePxMmStrategy
        for panel in detector:
            panel.set_px_mm_strategy(SimplePxMmStrategy())
            panel.set_mu(0)
            panel.set_thickness(0)
        assert detector[0].get_thickness() == 0

    alt_crystal = alt_expt.crystal

    beam = expt.beam
    spec = expt.imageset.get_spectrum(0)
    energies_raw = spec.get_energies_eV().as_numpy_array()
    weights_raw = spec.get_weights().as_numpy_array()
    from LS49.adse13_187.adse13_221.explore_spectrum import method3
    energies, weights, _ = method3(energies_raw, weights_raw,); weights = 5000000.*weights
    energies = list(energies); weights = list(weights)

    device_Id = 0
    if self.gpu_channels_singleton is not None:
      device_Id = self.gpu_channels_singleton.get_deviceID()

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
      mosaic_spread = params.mosaic_spread.value
      Ncells_abc = params.Nabc

      from LS49.adse13_187.cyto_batch import multipanel_sim
      JF16M_numpy_array, TIME_BG, TIME_BRAGG = multipanel_sim(
        CRYSTAL=alt_crystal, DETECTOR=detector, BEAM=beam,
        Famp = self.gpu_channels_singleton,
        energies=energies, fluxes=weights,
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

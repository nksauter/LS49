from __future__ import division, print_function
from time import time
"""
This is a specialization of cyto_batch.py script, but using the
exascale api implementation only.  The test is to run the simulation
with and without argchk and to assert that the results are the same.
"""
from LS49.adse13_187.cyto_batch import run_batch_job, multipanel_sim

import LS49.adse13_187.cyto_batch
def tst_one_monkeypatch(i_exp,spectra,Fmerge,gpu_channels_singleton,rank,params):
    print ("IN MONKEYPATCH")
    from simtbx.nanoBragg import utils
    from dxtbx.model.experiment_list import ExperimentListFactory
    import numpy as np

    print("Experiment %d" % i_exp, flush=True)

    outfile = "boop_%d.hdf5" % i_exp
    from LS49.adse13_187.case_data import retrieve_from_repo
    experiment_file = retrieve_from_repo(i_exp)
    cuda = True  # False  # whether to use cuda
    omp = False
    ngpu_on_node = 1 # 8  # number of available GPUs
    mosaic_spread = 0.07  # degrees
    mosaic_spread_samples = params.mosaic_spread_samples # number of mosaic blocks sampling mosaicity
    Ncells_abc = 30, 30, 10  # medians from best stage1
    ev_res = 1.5  # resolution of the downsample spectrum
    total_flux = 1e12  # total flux across channels
    beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
    spot_scale = 500. # 5.16324  # median from best stage1
    plot_spec = False  # plot the downsample spectra before simulating
    oversample = 1  # oversample factor, 1,2, or 3 probable enough
    panel_list = None  # integer list of panels, usefule for debugging
    rois_only = False  # only set True if you are running openMP, or CPU-only (i.e. not for GPU)
    include_background = params.include_background   # default is to add water background 100 mm thick
    verbose = 0  # leave as 0, unles debug
    flat = True  # enfore that the camera has 0 thickness
    #<><><><><><><><>
    # XXX new code
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

    if panel_list is None:
        panel_list = list(range(len(detector)))

    pids_for_rank = panel_list
    device_Id = 0
    if gpu_channels_singleton is not None:
      device_Id = gpu_channels_singleton.get_deviceID()

    print("Rank %d will use device %d" % (rank, device_Id))
    show_params = False
    time_panels = (rank == 0)

    mn_energy = (energies*weights).sum() / weights.sum()
    mn_wave = utils.ENERGY_CONV / mn_energy
    jf16m_numpy_array = {}
    from_gpu_amplitudes_cuda = {}
    print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
    print("\tBreakdown:")
    for shapetype in ["gauss_argchk","gauss"]:
      BEG=time()
      print (gpu_channels_singleton.get_deviceID(),"device",shapetype)
      Famp_is_uninitialized = ( gpu_channels_singleton.get_nchannels() == 0 )
      if Famp_is_uninitialized:
        F_P1 = Fmerge.expand_to_p1()
        for x in range(1):  # in this scenario, amplitudes are independent of lambda
          gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data())
      assert gpu_channels_singleton.get_nchannels() == 1

      JF16M_numpy_array, TIME_BG, TIME_BRAGG, _ = multipanel_sim(
        CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
        Famp = gpu_channels_singleton,
        energies=list(energies), fluxes=list(weights),
        background_wavelengths=[mn_wave], background_wavelength_weights=[1],
        background_total_flux=total_flux,background_sample_thick_mm=0.5,
        cuda=True,
        oversample=oversample, Ncells_abc=Ncells_abc,
        mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
        beamsize_mm=beamsize_mm,
        profile=shapetype,
        show_params=show_params,
        time_panels=time_panels, verbose=verbose,
        spot_scale_override=spot_scale,
        include_background=include_background)
      TIME_EXA = time()-BEG
      jf16m_numpy_array[shapetype]=JF16M_numpy_array
      from_gpu_amplitudes_cuda[shapetype]=TIME_BRAGG

      print("\t\tExascale: time for bkgrd sim: %.4fs; Bragg sim: %.4fs; total: %.4fs" % (TIME_BG, TIME_BRAGG, TIME_EXA))
    ratio = from_gpu_amplitudes_cuda["gauss"]/from_gpu_amplitudes_cuda["gauss_argchk"]
    print("<><><><><><><><><ratio<%.2f><><><><><><><><><><><><><><><><><><><>\n"%(
      ratio))

    # assertion on elapsed time:
    assert ratio > 1.0,"ratio is %.3f, experiment %d"%(ratio,i_exp)

    # assertion on equality:
    abs_diff = np.abs(jf16m_numpy_array["gauss"] - \
                      jf16m_numpy_array["gauss_argchk"]).max()
    assert np.allclose(jf16m_numpy_array["gauss"], \
                       jf16m_numpy_array["gauss_argchk"]), \
    "max per-pixel difference: %f photons, experiment %d"%(abs_diff,i_exp)

LS49.adse13_187.cyto_batch.tst_one = tst_one_monkeypatch

if __name__=="__main__":
  run_batch_job()
  print ("OK")

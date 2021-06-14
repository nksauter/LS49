from __future__ import division
import libtbx.load_env
import os
import pickle
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.array_family import flex
from simtbx.nanoBragg import utils
import numpy as np
import scipy

def tst_one_monkeypatch(i_exp,spectra,Fmerge,gpu_channels_singleton,rank,params):

  data_path = os.path.join(
      libtbx.env.find_in_repositories('ls49_big_data'),
      'adse13_228',
      'tst_whitelist'
  )
  data_root = 'idx-run_000452.JF07T32V01_master_00001'
  experiment_file = os.path.join(data_path, data_root+'_integrated.expt')
  refl_file = os.path.join(data_path, data_root+'_filtered.refl')
  spectrum_file = os.path.join(data_path, 'spectrum.pickle')
  mask_file = os.path.join(data_path, '4more.mask')
  results_file = os.path.join(data_path, 'results.npz')

  cuda = True  # False  # whether to use cuda
  mosaic_spread = 0.07  # degrees
  mosaic_spread_samples = params.mosaic_spread_samples # number of mosaic blocks sampling mosaicity
  Ncells_abc = 30, 30, 10  # medians from best stage1
  total_flux = 1e12  # total flux across channels
  beamsize_mm = 0.000886226925452758  # sqrt of beam focal area
  spot_scale = 500. # 5.16324  # median from best stage1
  oversample = 1  # oversample factor, 1,2, or 3 probable enough
  include_background = False
  verbose = 0  # leave as 0, unles debug
  shapetype = 'gauss_argchk'
  show_params=False
  time_panels = (rank==0)

  El = ExperimentListFactory.from_json_file(experiment_file, check_format=True)
  exper = El[0]
  crystal = exper.crystal
  detector = exper.detector
  beam = exper.beam

  with open(spectrum_file, 'rb') as f:
    energies = np.array(pickle.load(f))
    weights = np.array(pickle.load(f))
  mn_energy = (energies*weights).sum() / weights.sum()
  mn_wave = utils.ENERGY_CONV / mn_energy

  from LS49.adse13_187.adse13_221.basic_runner import basic_run_manager
  M = basic_run_manager.from_files(mask_file, refl_file, experiment_file)
  M.refl_analysis(experiment_file)
  M.get_trusted_and_refl_mask()

  # ignore Fmerge, make our own amplitudes
  from cctbx.crystal import symmetry
  from cctbx.miller import array as miller_array, set as miller_set
  uc = crystal.get_unit_cell()
  sg = crystal.get_space_group()
  MS = miller_set(
      symmetry(unit_cell=uc, space_group=sg), anomalous_flag=True,
      indices=M.refl_table['miller_index'].select(M.refl_table['spots_order'])
  )
  dummy_data = flex.double(
      [i%100 for i in range(len(M.refl_table['spots_order']))]
  )
  Fmerge_dummy = miller_array(MS, data=dummy_data)

  Famp_is_uninitialized = (gpu_channels_singleton.get_nchannels()==0)
  if Famp_is_uninitialized:
    F_P1 = Fmerge_dummy.expand_to_p1()
    for x in range(1):
      gpu_channels_singleton.structure_factors_to_GPU_direct_cuda(
          x, F_P1.indices(), F_P1.data()
      )

  mask_bool = M.monolithic_mask_whole_detector_as_1D_bool
  active_pixels = flex.int()
  for i, x in enumerate(mask_bool):
    if x: active_pixels.append(i)
  mask_int = active_pixels

  arr1, TIME_BG, TIME_BRAGG = multipanel_sim(
      CRYSTAL=crystal, DETECTOR=detector, BEAM=beam,
      Famp = gpu_channels_singleton,
      energies=list(energies), fluxes=list(weights),
      background_wavelengths=[mn_wave], background_wavelength_weights=[1],
      background_total_flux=total_flux,background_sample_thick_mm=0.5,
      cuda=cuda,
      oversample=oversample, Ncells_abc=Ncells_abc,
      mos_dom=mosaic_spread_samples, mos_spread=mosaic_spread,
      beamsize_mm=beamsize_mm,
      profile=shapetype,
      show_params=show_params,
      time_panels=time_panels, verbose=verbose,
      spot_scale_override=spot_scale,
      include_background=include_background,
      mask_file="")
  arr2, TIME_BG, TIME_BRAGG = multipanel_sim(
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
      include_background=include_background,
      mask_file=mask_int)
  arr1f = arr1.flatten()
  arr2f = arr2.flatten()
  dump_ref = False
  if dump_ref:
    arr2s = scipy.sparse.csc_matrix(arr2f)
    scipy.sparse.save_npz(results_file, arr2s)
    quit()
  ref = scipy.sparse.load_npz(results_file).toarray()[0]
  from libtbx.test_utils import approx_equal
  n_inside_mask = 0
  for i, val in enumerate(mask_bool):
    if val:
      assert approx_equal(arr2f[i], ref[i])
      assert approx_equal(arr1f[i], arr2f[i])
      n_inside_mask += 1
    elif i%97==0:
      assert approx_equal(arr2f[i], 0)
    if i%1000000==0:
      print('tested ', i, ' elements including ', n_inside_mask, ' inside mask')


from LS49.adse13_187.cyto_batch import run_batch_job, multipanel_sim
import LS49.adse13_187.cyto_batch
LS49.adse13_187.cyto_batch.tst_one = tst_one_monkeypatch

if __name__=='__main__':
  run_batch_job()

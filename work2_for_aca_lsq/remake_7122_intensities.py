from __future__ import print_function
from __future__ import division
from six.moves import cPickle as pickle
from six.moves import range

if __name__=="__main__":

  # %%% boilerplate specialize to packaged big data %%%
  import os
  from LS49.sim import step5_pad
  from LS49.sim import step4_pad
  from LS49.spectra import generate_spectra
  ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
  step5_pad.big_data = ls49_big_data
  step4_pad.big_data = ls49_big_data
  generate_spectra.big_data = ls49_big_data
  # %%%%%%

  from LS49.sim.util_fmodel import gen_fmodel
  from LS49.sim.step5_pad import data
  pdb_lines = data().get("pdb_lines")
  Fe_oxidized_model = data().get("Fe_oxidized_model")
  Fe_reduced_model = data().get("Fe_reduced_model")

  W2 = 12398.425/7122.

  GF = gen_fmodel(resolution=1.7,pdb_text=pdb_lines,algorithm="fft",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.params2.fmodel.b_sol = 46.
  GF.params2.structure_factors_accuracy.grid_resolution_factor = 1/5.
  GF.params2.mask.grid_step_factor = 10.
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_reduced = GF.get_intensities()
  # Einsle paper: Reduced form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(II) state (absorption at lower energy, reduced)

  W2i = W2_reduced.indices()
  with (open("confirm_intensities_flex_array.data","w")) as F:
    for iw in range(len(W2i)):
      print ("%20s, %10.2f"%(W2_reduced.indices()[iw],W2_reduced.data()[iw]), file=F)

  intensity_dict = {}
  for iw in range(len(W2i)):
    intensity_dict[W2_reduced.indices()[iw]] = W2_reduced.data()[iw]

  with (open("confirm_intensities_dict.pickle","wb")) as F:
    pickle.dump(intensity_dict, F, pickle.HIGHEST_PROTOCOL)

  with (open("confirm_C2_sfall_7122_amplitudes.pickle","wb")) as F:
    pickle.dump(GF.get_amplitudes(), F, pickle.HIGHEST_PROTOCOL)

  GF.make_P1_primitive()
  with (open("confirm_sfall_P1_7122_amplitudes.pickle","wb")) as F:
    pickle.dump(GF.get_amplitudes(), F, pickle.HIGHEST_PROTOCOL)

from __future__ import print_function
from __future__ import division
from six.moves import cPickle as pickle
from six.moves import range
from scitbx.array_family import flex

# %%% boilerplate specialize to packaged big data %%%
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%

from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")

def remake_intensities_at_energy(energy,FE1_model,FE2_model):
  from LS49.sim.util_fmodel import gen_fmodel

  W2 = 12398.425/float(energy)

  GF = gen_fmodel(resolution=1.7,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.params2.fmodel.b_sol = 46.
  GF.params2.structure_factors_accuracy.grid_resolution_factor = 1/5.
  GF.params2.mask.grid_step_factor = 10.
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=FE1_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=FE2_model,newvalue=W2)
  GF.make_P1_primitive()
  W2_reduced = GF.get_intensities()
  # Einsle paper: Reduced form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(II) state (absorption at lower energy, reduced)

  W2i = W2_reduced.indices()
  intensity_dict = {}
  for iw in range(len(W2i)):
    intensity_dict[W2_reduced.indices()[iw]] = W2_reduced.data()[iw]
  return intensity_dict

if __name__=="__main__":

  for filename,FE1,FE2 in [
    ("confirm_P1_range_reduced_intensities_dict.pickle", Fe_oxidized_model, Fe_reduced_model),
    ("confirm_P1_range_oxidized_intensities_dict.pickle", Fe_oxidized_model, Fe_oxidized_model),
    ("confirm_P1_range_metallic_intensities_dict.pickle", Fe_metallic_model, Fe_metallic_model),
    ("confirm_P1_range_swapreduced_intensities_dict.pickle", Fe_reduced_model, Fe_oxidized_model),
  ]:
    initial = remake_intensities_at_energy(7070.0,FE1,FE2)
    print ("%d keys in initial"%len(initial))

    result = {}
    for key in initial:
      result[key] = flex.double()
    print (filename)
    for incr in range(100):
      energy = 7070.5 + incr
      more = remake_intensities_at_energy(energy,FE1,FE2)
      print (energy, "with %d keys"%(len(more)))
      for key in more:
        if key in initial:
          result[key].append(more[key])
    exit()
    with (open(filename,"wb")) as F:
      pickle.dump(result, F, pickle.HIGHEST_PROTOCOL)


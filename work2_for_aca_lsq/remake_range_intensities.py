from __future__ import print_function
from __future__ import division
from six.moves import cPickle as pickle
from six.moves import range
from scitbx.array_family import flex

def remake_intensities_at_energy(energy):
  from LS49.sim.util_fmodel import gen_fmodel
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model
  from LS49.sim.fdp_plot import george_sherrell
  #Fe_oxidized_model = george_sherrell("data_sherrell/pf-rd-ox_fftkk.out")
  #Fe_reduced_model = george_sherrell("data_sherrell/pf-rd-red_fftkk.out")
  Fe_metallic_model = george_sherrell("data_sherrell/Fe_fake.dat")

  W2 = 12398.425/float(energy)

  GF = gen_fmodel(resolution=1.7,pdb_text=pdb_lines,algorithm="fft",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.params2.fmodel.b_sol = 46.
  GF.params2.structure_factors_accuracy.grid_resolution_factor = 1/5.
  GF.params2.mask.grid_step_factor = 10.
  GF.reset_wavelength(W2)
  #GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2) # Einsle ox and reduced
  #GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_oxidized_model,newvalue=W2) # Einsle oxidized
  #GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2) # Einsle reduced
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_metallic_model,newvalue=W2) # Einsle metallic
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_metallic_model,newvalue=W2) # Einsle metallic
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

  initial = remake_intensities_at_energy(7070.0)
  print ("%d keys in initial"%len(initial))

  result = {}
  for key in initial:
    result[key] = flex.double()

  for incr in range(100):
    energy = 7070.5 + incr
    more = remake_intensities_at_energy(energy)
    print (energy, "with %d keys"%(len(more)))
    for key in more:
      if key in initial:
        result[key].append(more[key])

  #with (open("confirm_P1_range_intensities_dict.pickle","wb")) as F: # Einsle reduced
  #with (open("confirm_P1_range_oxidized_intensities_dict.pickle","wb")) as F: # Einsle oxidized
  with (open("confirm_P1_range_metallic_intensities_dict.pickle","wb")) as F: # Einsle metallic
    pickle.dump(result, F, pickle.HIGHEST_PROTOCOL)

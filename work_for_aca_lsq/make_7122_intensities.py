from __future__ import print_function
from __future__ import division
from cctbx.array_family import flex
import glob
import scitbx
from six.moves import cPickle as pickle
from six.moves import range

if __name__=="__main__":

  from LS49.sim.util_fmodel import gen_fmodel
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model

  W2 = 12398.425/7122.

  GF = gen_fmodel(resolution=1.9,pdb_text=pdb_lines,algorithm="fft",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_reduced = GF.get_intensities()
  # Einsle paper: Reduced form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(II) state (absorption at lower energy, reduced)

  from cctbx import miller
  W2i = W2_reduced.indices()
  with (open("debug26.data","w")) as F:
    for iw in range(len(W2i)):
      print ("%20s, %10.2f"%(W2_reduced.indices()[iw],W2_reduced.data()[iw]), file=F)

  intensity_dict = {}
  for iw in range(len(W2i)):
    intensity_dict[W2_reduced.indices()[iw]] = W2_reduced.data()[iw]

  with (open("debug26_intensities.pickle","wb")) as F:
    pickle.dump(intensity_dict, F, pickle.HIGHEST_PROTOCOL)

  with (open("sfall_7122_amplitudes.pickle","wb")) as F:
    pickle.dump(GF.get_amplitudes(), F, pickle.HIGHEST_PROTOCOL)

  GF.make_P1_primitive()
  with (open("sfall_P1_7122_amplitudes.pickle","wb")) as F:
    pickle.dump(GF.get_amplitudes(), F, pickle.HIGHEST_PROTOCOL)

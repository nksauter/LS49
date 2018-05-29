from __future__ import print_function
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
  # Einsle paper: Reduced form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(II) state (absorption at lower energy, reduced)

  per_HKL_I = {}
  for iE,Energy in enumerate(xrange(7090,7151)):
    if Energy%10==0: print (Energy)
    W = 12398.425/Energy
    GF.reset_wavelength(W)
    GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W)
    GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W)
    this_energy_Fmodel = GF.get_intensities()
    millers = this_energy_Fmodel.indices()
    data = this_energy_Fmodel.data()
    if iE==0:
      for key in millers: per_HKL_I[key]=flex.double()
    for ikey in xrange(len(millers)):
      if millers[ikey] in per_HKL_I:
        per_HKL_I[millers[ikey]].append(data[ikey])

  with (open("debug26_range_intensities.pickle","wb")) as F:
    pickle.dump(per_HKL_I, F, pickle.HIGHEST_PROTOCOL)

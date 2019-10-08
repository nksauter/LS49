from __future__ import division, print_function
from six.moves import cPickle
import os
import six

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment

def get_pdb_lines():
  return open(os.path.join(ls49_big_data,"1m2a.pdb"),"r").read()

def single_wavelength_fmodel(create):
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_images(20,energy=7120.,total_flux=1e12)
  wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength

  direct_algo_res_limit = 1.7
  from LS49.sim.util_fmodel import gen_fmodel
  for flag in [True,False]:
    GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=get_pdb_lines(),algorithm="fft",wavelength=wavelength_A)
    GF.set_k_sol(0.435)
    if flag:  GF.make_P1_primitive()
    sfall_main = GF.get_amplitudes()
    sfall_main.show_summary(prefix = "Amplitudes used ")
    if create: # write the reference for the first time
      cPickle.dump(sfall_main,
      open(os.path.join(ls49_big_data,"reference","sf_reference_cb_to_P1_%s"%(str(flag))),"wb"),cPickle.HIGHEST_PROTOCOL)
    else: # read the reference and assert sameness to sfall_main

      if six.PY3:
        sfall_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference","sf_reference_cb_to_P1_%s"%(str(flag))),"rb"), encoding="bytes")
        from LS49.tests.tst_sf_energies import fix_unpickled_attributes
        fix_unpickled_attributes(sfall_ref)
      else:
        sfall_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference","sf_reference_cb_to_P1_%s"%(str(flag))),"rb"))



      T = sfall_main; S = sfall_ref
      assert S.space_group() == T.space_group()
      assert S.unit_cell().parameters() == T.unit_cell().parameters()
      assert S.indices() == T.indices(); assert S.data() == T.data()
    print()

# ideally, do a calculation to prove that the duplicated structure factors in the
# P1 cell are equivalent to those in C2, and / or try to attribute the difference
# to the bulk solvent model rather then the site scatterer model.

if __name__=="__main__":
  single_wavelength_fmodel(create=False)
  print("OK")

from __future__ import division, print_function
from six.moves import cPickle
import os
from LS49.spectra import generate_spectra

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
def get_results():
  return cPickle.load(open(os.path.join(ls49_big_data,"data/spectra209.pickle"),"rb"))
generate_spectra.get_results = get_results

def get_pdb_lines():
  return open(os.path.join(ls49_big_data,"1m2a.pdb"),"r").read()

ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment

from LS49.sim.fdp_plot import george_sherrell
Fe_oxidized_model = george_sherrell(os.path.join(ls49_big_data,"data_sherrell/pf-rd-ox_fftkk.out"))
Fe_reduced_model = george_sherrell(os.path.join(ls49_big_data,"data_sherrell/pf-rd-red_fftkk.out"))

def channel_wavelength_fmodel(create):
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  iterator = SS.generate_recast_renormalized_images(20,energy=7120.,total_flux=1e12)
  wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength

  direct_algo_res_limit = 1.7
  from LS49.sim.util_fmodel import gen_fmodel

  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=get_pdb_lines(),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()
  for x in range(10,len(flux),5):
    print("+++++++++++++++++++++++++++++++++++++++ Wavelength",x)
    GF.reset_wavelength(wavelength_A)
    GF.reset_specific_at_wavelength(
                   label_has="FE1",tables=Fe_oxidized_model,newvalue=wavelength_A)
    GF.reset_specific_at_wavelength(
                   label_has="FE2",tables=Fe_reduced_model,newvalue=wavelength_A)
    sfall_channel = GF.get_amplitudes()
    filename = "sf_reference_channel_%s"%("%03d"%x)
    if create: # write the reference for the first time
      cPickle.dump(sfall_channel,
      open(os.path.join(ls49_big_data,"reference",filename),"wb"),cPickle.HIGHEST_PROTOCOL)
    else: # read the reference and assert sameness to sfall_channel
      sfall_ref = cPickle.load(open(os.path.join(ls49_big_data,"reference",filename),"rb"))
      T = sfall_channel; S = sfall_ref
      assert S.space_group() == T.space_group()
      assert S.unit_cell().parameters() == T.unit_cell().parameters()
      assert S.indices() == T.indices(); assert S.data() == T.data()

if __name__=="__main__":
  channel_wavelength_fmodel(create=False)
  print("OK")

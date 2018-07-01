from __future__ import print_function
from __future__ import division
from cctbx.array_family import flex
from six.moves import cPickle as pickle
from six.moves import range
from LS49.work2_for_aca_lsq.make_model_obs_28 import GM # implicit import
from matplotlib import pyplot as plt


def test_Gi_factor(G):
  from LS49.spectra.generate_spectra import spectra_simulation
  SS = spectra_simulation()
  x = flex.double()
  y = flex.double()
  for key in G.images_strong:
    print (key,"of",len(G.images_strong),G.images_Gi[key])
    #trying here to plot the Gi against the integrated spectrum for each event.
    iterator = SS.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)
    wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength
    total_flux = flex.sum(flux)
    x.append(total_flux)
    y.append(G.images_Gi[key])
  from matplotlib import pyplot as plt
  plt.plot(x,y,"r.")
  plt.show()


def get_C2_pdb_structure(resolution,wavelength):
    from LS49.sim.util_fmodel import gen_fmodel
    pdb_lines = open("./1m2a.pdb","r").read()
    return gen_fmodel(
         resolution=resolution, pdb_text=pdb_lines,
         wavelength=wavelength, algorithm="fft"
    )

def reproduce_C2_model(resolution,**kwargs):
  from LS49.sim.step5_pad import Fe_oxidized_model,Fe_reduced_model
  W2 = 12398.425/7122.
  GF = get_C2_pdb_structure(resolution,wavelength=W2)
  GF.params2.fmodel.b_sol = 46.
  GF.set_k_sol(0.435)
  # take defaults: b_sol=46., grid_resolution_factor=1/3., grid_step_factor=4.
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  sfall = GF.get_amplitudes()
  sfallf = sfall.data(); sfallidx = sfall.indices()
  intensities = {}
  for N in range(len(sfallidx)):
    intensities[sfallidx[N]]=sfallf[N] * sfallf[N]
  return intensities

if __name__=="__main__":

  # pickle contains the results of the "partiality" simulation (top-hat incident light, debug26_intensities)
  # actually, no, this particular simulation (using current code in util_partiality.py) apparently
  # uses sfall_P1_7122_amplitudes.pickle then of course inserts this intensity into the dictionary key "simtbx_intensity".
  # by the way, here are the parameters for sfall_P1_7122_amplitudes
  #   W2 = 12398.425/7122.
  #   GF = gen_fmodel(resolution=1.9,pdb_text=pdb_lines,algorithm="fft",wavelength=W2)
  #   GF.set_k_sol(0.435)
  #   GF.reset_wavelength(W2)
  #   GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  #   GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  # make_model_obs_28.py scores simulation against experiment, selects reflections with high CC,
  # then writes them to this single pickle file.
  with open("make_model_obs_28.pickle","rb") as VVV:
    G = pickle.load(VVV)

  #with (open("debug26_range_intensities.pickle","rb")) as F:
  with (open("confirm_range_intensities_dict.pickle","rb")) as F:
    per_HKL_I = pickle.load(F)

  #with (open("debug26_intensities.pickle","rb")) as F:
  #  per_HKL_I_7122 = pickle.load(F)
  # let's use the experience of C2_to_P1_demo.py to work up a simulation of the "simtbx_intensity" really
  # used for the may27_gen_data_mpi.py partiality simulation
  per_HKL_I_7122 = reproduce_C2_model(resolution = 1.7)

  def method2_include_pehHKL_I_explicitly():
   # There's no reason why we can't get the Gi's by analytical least squares
   skeys = list(G.images_strong.keys())
   skeys.sort()
   for key in skeys:
    print ("image",key)
    numerator = 0.; denominator = 0.
    nkeys = len(G.images_strong[key])

    for ikey, HKL in enumerate(G.images_strong[key]):
      MD = G.images_strong[key][HKL]
      #from IPython import embed; embed()

      terms1 = MD["model"] * per_HKL_I[HKL] / MD["simtbx_intensity"]
      terms2 = terms1 * terms1
      terms0 = MD["obs"] * terms1
      numerator+=flex.sum(terms0)
      denominator+=flex.sum(terms2)
    G.images_Gi[key]=numerator/denominator

    for ikey, HKL in enumerate(G.images_strong[key]):
      plt.subplot(nkeys,1,1+ikey)
      MD = G.images_strong[key][HKL]
      assert len(MD["obs"])==61
      print (HKL,MD,"7122 lookup",per_HKL_I_7122[HKL],per_HKL_I_7122[HKL]/MD["simtbx_intensity"])
      plt.plot(range(7090,7151),(G.images_Gi[key]) * MD["model"] * per_HKL_I[HKL] / MD["simtbx_intensity"],"b-")
      plt.plot(range(7090,7151),MD["obs"],"r-")
#start here.  can we shwo we are actually at a minimum of teh target function, considering it is LSQ?
    plt.show()
    if key%100==0: print (key, "Gi:", G.images_Gi[key])

  method2_include_pehHKL_I_explicitly()
  print ("FINISHED LSQ GI, now testing against total flux")
  test_Gi_factor(G)
  """start here.  why are our Gi factors so screwed up?"""

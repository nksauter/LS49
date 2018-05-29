from __future__ import print_function
from cctbx.array_family import flex
import glob
import scitbx
from six.moves import cPickle as pickle
from six.moves import range
from seriously_deal_with_f_derivatives import eV_as_angstroms,at_one_eV
from make_model_obs import GM
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
    wavlen, flux, wavelength_A = iterator.next() # list of lambdas, list of fluxes, average wavelength
    total_flux = flex.sum(flux)
    x.append(total_flux)
    y.append(G.images_Gi[key])
  from matplotlib import pyplot as plt
  plt.plot(x,y,"r.")
  plt.show()

if __name__=="__main__":

  with open("make_model_obs_28.pickle","rb") as VVV:
    G = pickle.load(VVV)

  with (open("debug26_range_intensities.pickle","rb")) as F:
    per_HKL_I = pickle.load(F)

  with (open("debug26_intensities.pickle","rb")) as F:
    per_HKL_I_7122 = pickle.load(F)

  def method1_include_pehHKL_I_explicitly():
   # There's no reason why we can't get the Gi's by analytical least squares
   for key in G.images_strong:
    print (key, G.images_strong[key])
    numerator = 0.; denominator = 0.
    nkeys = len(G.images_strong[key])

    for ikey, HKL in enumerate(G.images_strong[key]):
      plt.subplot(nkeys,1,1+ikey)
      MD = G.images_strong[key][HKL]
      plt.plot(range(7090,7151),per_HKL_I[HKL]*MD["model"],"k-")
      plt.plot(range(7090,7151),1E14*MD["obs"],"r-")
      terms1 = G.images_strong[key][HKL]["model"] * per_HKL_I[HKL]
      terms2 = terms1 * terms1
      terms0 = G.images_strong[key][HKL]["obs"] * terms1
      numerator+=flex.sum(terms0)
      denominator+=flex.sum(terms2)
    plt.show()
    G.images_Gi[key]=numerator/denominator
    if key%100==0: print (key, "Gi:", G.images_Gi[key])

  def method2_include_pehHKL_I_explicitly():
   # There's no reason why we can't get the Gi's by analytical least squares
   for key in G.images_strong:
    numerator = 0.; denominator = 0.
    nkeys = len(G.images_strong[key])

    for ikey, HKL in enumerate(G.images_strong[key]):
      plt.subplot(nkeys,1,1+ikey)
      MD = G.images_strong[key][HKL]
      plt.plot(range(7090,7151),MD["model"] * per_HKL_I[HKL] / per_HKL_I_7122[HKL],"k-")
      plt.plot(range(7090,7151),1E10*MD["obs"],"r-")
      terms1 = MD["model"] * per_HKL_I[HKL] / per_HKL_I_7122[HKL]
      terms2 = terms1 * terms1
      terms0 = MD["obs"] * terms1
      numerator+=flex.sum(terms0)
      denominator+=flex.sum(terms2)
    plt.show()
    G.images_Gi[key]=numerator/denominator
    if key%100==0: print (key, "Gi:", G.images_Gi[key])

  method2_include_pehHKL_I_explicitly()
  print ("FINISHED LSQ GI, now testing against total flux")
  test_Gi_factor(G)
  """start here.  why are our Gi factors so screwed up?"""


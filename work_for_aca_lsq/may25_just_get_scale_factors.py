from __future__ import print_function
from __future__ import division
from cctbx.array_family import flex
import glob
from six.moves import cPickle as pickle
from six.moves import range
"""Focus on one energy at a time.  Assume we already know the image scale factors.
Only thing left to do is optimize the four parameters FE1 f'f" and FE2 f'f"
"""
class GM (object):
 def generate_millers(self):
  #globv = "../work/dataX0[0-1].pickle" #test
  globv = "../work/dataX*.pickle" # production
  self.asu = {}
  self.icount = 0
  self.images_all = 0
  self.images_strong = {}
  self.images_Gi = {}
  for filen in glob.glob(globv):

   with open(filen,"rb") as V:
    while 1:
     try:
      image = pickle.load(V)

      self.images_all+=1
      highcc = flex.double(image["cc"]) > 0.80 #was 0.7 and minimum 4
      if highcc.count(True)<3: continue
      imdict = dict()
      self.images_strong[image["image"]]=imdict
      for i in range(len(image["cc"])):
        if image["cc"][i]<0.8: continue
        self.icount+=1
        imdict[image["millers"][i]]=dict(model=image["model"][i], obs=flex.double(image["obs"][i]))
        print (filen,self.icount,"CC>70%%: %20s %5.2f"%(image["millers"][i],image["cc"][i]))
        #from matplotlib import pyplot as plt
        #model = image["model"][i]
        #obs = image["obs"][i]
        #plt.plot(range(len(model)),model,"k-")
        #plt.plot(range(len(obs)),1.E10*flex.double(obs),"r-")
        #plt.show()
        self.asu[image["millers"][i]]=1

        M = image["millers"][i]
        foundM = False
        with open("debug26.data","r") as VVV:
          for line in VVV:
            SS = str(M)
            if SS in line: foundM = True; break
        if foundM:  print ("found",M,line)
        else: print ("NOT F",M)

        yield image["millers"][i]

     except EOFError:
      break
from may25_ratio_approach import lbfgs_fpfdp_fit

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
  #from matplotlib import pyplot as plt
  #plt.plot(x,y,"r.")
  #plt.show()

if __name__=="__main__":
  M = flex.miller_index()
  G = GM()
  for i in G.generate_millers():
    M.append(i)
  print ("%d Bragg spots measured"%len(M))

  print ("%d Unique Miller indices"%(len(G.asu)))
  exit()
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
  matches = miller.match_indices(M,W2i)
  sel0 = flex.size_t([p[0] for p in matches.pairs()])
  sel1 = flex.size_t([p[1] for p in matches.pairs()])
  print ("matches",len(sel0))
  #sel0unique is a dictionary keyed by unique Miller indices.
  sel0unique={}
  for item in sel1:
    sel0unique[W2i[item]]=item
  print ("unique",len(sel0unique))
  print ("total images %d, strong %d"%(G.images_all,len(G.images_strong)))
  exit("done outputting intensities")
  #doublecheck this
  #for si in xrange(len(sel0)):
    #print (si, M[sel0[si]], W2i[sel1[si]], W2i[sel0unique[W2i[sel1[si]]]])

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
  exit()
  # There's no reason why we can't get the Gi's by analytical least squares
  for key in G.images_strong:
    #print (key, G.images_strong[key])
    numerator = 0.; denominator = 0.
    for HKL in G.images_strong[key]:
      #from IPython import embed; embed()
      terms1 = G.images_strong[key][HKL]["model"] * per_HKL_I[HKL]
      terms2 = terms1 * terms1
      terms0 = G.images_strong[key][HKL]["obs"] * terms1
      numerator+=flex.sum(terms0)
      denominator+=flex.sum(terms2)
    G.images_Gi[key]=numerator/denominator
    if key%100==0: print (key, "Gi:", G.images_Gi[key])

  test_Gi_factor(G)
  """start here.  why are our Gi factors so screwed up?"""
  exit()

  print("Got Gi scale factors")
  from IPython import embed; embed()
  # Now loop over energies and try to focus on the fpfdp problem
  result_energies = flex.double()
  result_FE1_fpfdp = flex.vec2_double()
  result_FE2_fpfdp = flex.vec2_double()
  for iE,Energy in enumerate(xrange(7110,7131)):
    Eidx = iE+20 # index into the G arrays, for that particular energy
    try:
      A = lbfgs_fpfdp_fit(energy = Energy, all_data = G, Eidx = Eidx)
      result_energies.append(Energy)
      result_FE1_fpfdp.append((A.a[0],A.a[1]))
      result_FE2_fpfdp.append((A.a[2],A.a[3]))
    except Exception:
      pass

  from matplotlib import pyplot as plt
  from LS49.sim.fdp_plot import george_sherrell
  GS = george_sherrell("data_sherrell/pf-rd-ox_fftkk.out")
  GS.plot_them(plt,f1="b.",f2="b.")
  GS.plot_them(plt,f1="b-",f2="b-")
  GS = george_sherrell("data_sherrell/pf-rd-red_fftkk.out")
  GS.plot_them(plt,f1="r.",f2="r.")
  GS.plot_them(plt,f1="r-",f2="r-")
  plt.plot(result_energies, result_FE1_fpfdp.parts()[0], "c-")
  plt.plot(result_energies, result_FE1_fpfdp.parts()[1], "c-")
  plt.plot(result_energies, result_FE2_fpfdp.parts()[0], "m-")
  plt.plot(result_energies, result_FE2_fpfdp.parts()[1], "m-")
  plt.xlabel('Energy (eV)')
  plt.xlim([7088,7152])
  plt.ylim([-8.2,4.2])
  plt.show()
  exit("STOP")

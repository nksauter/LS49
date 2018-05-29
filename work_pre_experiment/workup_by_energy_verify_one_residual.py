from __future__ import print_function
from cctbx.array_family import flex
import pickle,glob
import scitbx
from seriously_deal_with_f_derivatives import eV_as_angstroms,at_one_eV
"""Focus on one energy at a time.  Assume we already know the image scale factors.
Only thing left to do is optimize the four parameters FE1 f'f" and FE2 f'f"
"""
class GM (object):
 def generate_millers(self):
  globv = "dataX0[0-1].pickle"
  #globv = "dataX*.pickle"
  self.asu = {}
  self.icount = 0
  self.images_all = 0
  self.images_strong = {}
  self.images_Gi = {}
  for filen in glob.glob(globv):

    V = open(filen,"rb")
    while 1:
     try:
      image = pickle.load(V)

      self.images_all+=1
      highcc = flex.double(image["cc"]) > 0.70
      if highcc.count(True)<4: continue
      imdict = dict()
      self.images_strong[image["image"]]=imdict
      for i in range(len(image["cc"])):
        if image["cc"][i]<0.7: continue
        self.icount+=1
        imdict[image["millers"][i]]=dict(model=image["model"][i], obs=flex.double(image["obs"][i]))
        print (filen,self.icount,"CC>70%%: %20s %5.2f"%(image["millers"][i],image["cc"][i]))
        self.asu[image["millers"][i]]=1
        yield image["millers"][i]

     except EOFError,e:
      break
def lsq_target_function(title,label_table,images_Gi,genfmodel,genmiller):

  per_energy_I = {}
  per_HKL_I = {}
  for iE,Energy in enumerate(xrange(7090,7151)):
    if Energy%10==0: print (Energy)
    W = 12398.425/Energy
    genfmodel.reset_wavelength(W)
    for sitekey in label_table:
      genfmodel.reset_specific_at_wavelength(label_has=sitekey,tables=label_table[sitekey],newvalue=W)
    W_state = genfmodel.get_intensities()
    W2i = W_state.indices()
    matches = miller.match_indices(genmiller.asu.keys(),W2i)
    sel0 = flex.size_t([p[0] for p in matches.pairs()])
    sel1 = flex.size_t([p[1] for p in matches.pairs()])
    selected_I = W_state.data().select(sel1)
    selected_H = W_state.indices().select(sel1)
    if iE==0:
      for key in selected_H: per_HKL_I[key]=flex.double()
    for ikey in xrange(len(selected_H)):
      per_HKL_I[selected_H[ikey]].append(selected_I[ikey])
    per_energy_I[Energy] = selected_I
  # that gives all intensities at all energies

  # Evaluate the overall target function for the reduced state
  LSQ = 0.
  for key in G.images_strong:
    for HKL in G.images_strong[key]:
      terms1 = G.images_strong[key][HKL]["model"] * per_HKL_I[HKL]
      termsm = images_Gi[key] * terms1
      termsr = G.images_strong[key][HKL]["obs"] - termsm
      sumrsq = flex.sum(termsr*termsr)
      LSQ += sumrsq
  print ("The residual for %s is "%title, LSQ)

class lbfgs_fpfdp_fit:
  def __init__(self,energy,all_data,Eidx):
    print ("LBFGS energy",energy)
    self.G = all_data
    self.energy = energy
    self.Eidx = Eidx #index for this energy into the arrays of self.G
    self.n = 4
    self.x = flex.double([-7.,2.,-7.,2.]) # for now, always initialize to same
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]")

  def lsq_target_function(self,title,debug):

    per_HKL_I = {}
    miller_lookup = {}

    for iH, asuHKL in enumerate(self.miller):
      per_HKL_I[asuHKL] = self.intensity[iH]
      miller_lookup[asuHKL] = iH

    # Evaluate the overall target function for the reduced state
    LSQ = 0.
    self.debug_residuals = flex.double()
    self.debug_dresiduals = flex.double()
    self.gvec = flex.double(4)
    for key in self.G.images_strong:  # 733 images
      for HKL in self.G.images_strong[key]:  # should be 3414 terms visited in the double loop
        terms1 = self.G.images_strong[key][HKL]["model"][self.Eidx] * per_HKL_I[HKL]
        termsm = self.G.images_Gi[key] * terms1
        termsr = self.G.images_strong[key][HKL]["obs"][self.Eidx] - termsm
        self.debug_residuals.append(termsr)
        sumrsq = termsr*termsr
        LSQ += 0.5 * sumrsq
        gvec_fac = -2.0 * termsr * self.G.images_Gi[key] * self.G.images_strong[key][HKL]["model"][self.Eidx]
        termsdr = -2.0 * self.G.images_Gi[key] * self.G.images_strong[key][HKL]["model"][self.Eidx]
        self.gvec[0] += gvec_fac * self.d1[miller_lookup[HKL]]
        self.gvec[1] += gvec_fac * self.d2[miller_lookup[HKL]]
        self.gvec[2] += gvec_fac * self.d3[miller_lookup[HKL]]
        self.gvec[3] += gvec_fac * self.d4[miller_lookup[HKL]]
        self.debug_dresiduals.append(termsdr * self.d2[miller_lookup[HKL]])
    return LSQ

  def functional(self,values):
    self.miller, \
    self.intensity, \
    self.d1, self.d2, self.d3, self.d4 = at_one_eV(eV = self.energy, values = values)


    f = self.lsq_target_function(title="Internal form",debug="Residual")
    residuals0 = self.debug_residuals
    dr0 = self.debug_dresiduals

    v2 = [values[0]+0.0000, values[1]+0.0001, values[2], values[3]]
    self.miller, \
    self.intensity, \
    d1, d2, d3, d4 = at_one_eV(eV = self.energy, values = flex.double(v2))

    ff = self.lsq_target_function(title="Internal form",debug="Intensity")
    residuals1 = self.debug_residuals

    inn = len(residuals0)
    A = (residuals1-residuals0)/0.0001
    for ix in xrange(inn):
      print (ix,A[ix],dr0[ix])
    #Verified derivatives for each residual
start here:
    # Now what about each term in the functional?

    return f

  def gvec_callable(self,values):
    return self.gvec

  def compute_functional_and_gradients(self):
    self.a = self.x
    f = self.functional(self.x)
    self.print_step("LBFGS stp",f)
    g = self.gvec_callable(self.x)
    return f, g

if __name__=="__main__":
  M = flex.miller_index()
  G = GM()
  for i in G.generate_millers():
    M.append(i)
  print ("%d Bragg spots measured"%len(M))

  print ("%d Unique Miller indices"%(len(G.asu)))

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

  #doublecheck this
  #for si in xrange(len(sel0)):
    #print (si, M[sel0[si]], W2i[sel1[si]], W2i[sel0unique[W2i[sel1[si]]]])

  per_energy_I = {}
  per_HKL_I = {}
  for iE,Energy in enumerate(xrange(7090,7151)):
    if Energy%10==0: print (Energy)
    W = 12398.425/Energy
    GF.reset_wavelength(W)
    GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W)
    GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W)
    W_reduced = GF.get_intensities()
    W2i = W_reduced.indices()
    matches = miller.match_indices(G.asu.keys(),W2i)
    sel0 = flex.size_t([p[0] for p in matches.pairs()])
    sel1 = flex.size_t([p[1] for p in matches.pairs()])
    selected_I = W_reduced.data().select(sel1)
    selected_H = W_reduced.indices().select(sel1)
    if iE==0:
      for key in selected_H: per_HKL_I[key]=flex.double()
    for ikey in xrange(len(selected_H)):
      per_HKL_I[selected_H[ikey]].append(selected_I[ikey])
    per_energy_I[Energy] = selected_I
  # that gives all intensities at all energies

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

  print("Got Gi scale factors")
  # Now loop over energies and try to focus on the fpfdp problem
  result_energies = flex.double()
  result_FE1_fpfdp = flex.vec2_double()
  result_FE2_fpfdp = flex.vec2_double()
  for iE,Energy in enumerate(xrange(7110,7131)):
    if Energy !=7122 and Energy !=7123: continue
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

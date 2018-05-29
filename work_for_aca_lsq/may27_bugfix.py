from __future__ import print_function
from cctbx.array_family import flex
import pickle,glob
import scitbx.lbfgs
from seriously_deal_with_f_derivatives import eV_as_angstroms,at_one_eV
"""Focus on one energy at a time.  Assume we already know the image scale factors.
Only thing left to do is optimize the four parameters FE1 f'f" and FE2 f'f"
"""
from make_model_obs_28 import GM

class lbfgs_fpfdp_fit:
  def __init__(self,energy,all_data,Eidx):
    print ("LBFGS energy",energy, "with index",Eidx)
    self.G = all_data
    self.energy = energy
    self.Eidx = Eidx #index for this energy into the arrays of self.G
    self.n = 4
    self.x = flex.double([-7.,2.,-7.,2.]) # for now, always initialize to same
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=False,
        drop_convergence_test_max_drop_eps=1.E-1)
    )
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%10.4f"%a for a in self.x]),"]",", %d refl,%d images"%(self.nvisited,len(self.G.images_strong)))

  def lsq_target_function(self,title,debug=None):
     
    per_step_per_HKL_I = {}
    miller_lookup = {}

    for iH, asuHKL in enumerate(self.miller):
      per_step_per_HKL_I[asuHKL] = self.intensity[iH]
      miller_lookup[asuHKL] = iH

    # Evaluate the overall target function for the reduced state
    LSQ = 0.
    self.debug_terms = flex.double()
    self.debug_dterms = flex.double()
    self.gvec = flex.double(4)
    self.nvisited = 0
    for key in self.G.images_strong:  # 733 images
      for HKL in self.G.images_strong[key]:  # should be 3414 terms visited in the double loop
        self.nvisited+=1
        terms1 = self.G.images_strong[key][HKL]["model"][self.Eidx] * per_step_per_HKL_I[HKL] / per_HKL_I_7122[HKL]
        termsm = self.G.images_Gi[key] * terms1
        termsr = self.G.images_strong[key][HKL]["obs"][self.Eidx] - termsm
        sumrsq = termsr*termsr
        if debug is not None: self.debug_terms.append(0.5*sumrsq)
        LSQ += 0.5 * sumrsq
        gvec_fac = -2.0 * termsr * self.G.images_Gi[key] * self.G.images_strong[key][HKL]["model"][self.Eidx] / per_HKL_I_7122[HKL]
        if debug is not None: termsdr = -2.0 * self.G.images_Gi[key] * self.G.images_strong[key][HKL]["model"][self.Eidx]
        self.gvec[0] += gvec_fac * self.d1[miller_lookup[HKL]]
        self.gvec[1] += gvec_fac * self.d2[miller_lookup[HKL]]
        self.gvec[2] += gvec_fac * self.d3[miller_lookup[HKL]]
        self.gvec[3] += gvec_fac * self.d4[miller_lookup[HKL]]
        if debug is not None: self.debug_dterms.append(gvec_fac * self.d1[miller_lookup[HKL]])
    return LSQ

  def functional(self,values):
    self.miller, \
    self.intensity, \
    self.d1, self.d2, self.d3, self.d4 = at_one_eV(eV = self.energy, values = values)


    f = self.lsq_target_function(title="Internal form")
    return f
    terms0 = self.debug_terms
    dt0 = self.debug_dterms

    v2 = [values[0]+0.0001, values[1]+0.0000, values[2], values[3]]
    self.miller, \
    self.intensity, \
    d1, d2, d3, d4 = at_one_eV(eV = self.energy, values = flex.double(v2))

    ff = self.lsq_target_function(title="Internal form",debug="Intensity")
    terms1 = self.debug_terms

    inn = len(terms0)
    A = (terms1-terms0)/0.0001
    for ix in xrange(inn):
      print (ix,A[ix],dt0[ix]) 
    # Verified each term in the functional?

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
  with open("make_model_obs_28.pickle","rb") as VVV:
    G = pickle.load(VVV)

  with (open("debug26_range_intensities.pickle","rb")) as F:
    per_HKL_I = pickle.load(F)

  with (open("debug26_intensities.pickle","rb")) as F:
    per_HKL_I_7122 = pickle.load(F)

  # There's no reason why we can't get the Gi's by analytical least squares
  for key in G.images_strong:
    numerator = 0.; denominator = 0.
    nkeys = len(G.images_strong[key])

    for ikey, HKL in enumerate(G.images_strong[key]):
      MD = G.images_strong[key][HKL]
      terms1 = G.images_strong[key][HKL]["model"] * per_HKL_I[HKL] / per_HKL_I_7122[HKL]
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
    Eidx = iE+20 # index into the G arrays, for that particular energy
    if 1:#try:
      A = lbfgs_fpfdp_fit(energy = Energy, all_data = G, Eidx = Eidx)
      result_energies.append(Energy)
      result_FE1_fpfdp.append((A.a[0],A.a[1]))
      result_FE2_fpfdp.append((A.a[2],A.a[3]))
    if 0:#except Exception:
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


  


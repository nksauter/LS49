from __future__ import print_function, absolute_import
from __future__ import division
from six.moves import range
from cctbx.array_family import flex
import pickle,glob

class GM (object):
 def generate_millers(self):
  globv = "dataX0[0-7].pickle"
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
      print (image["image"])

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

     except EOFError as e:
      break
def lsq_target_function(title,label_table,images_Gi,genfmodel,genmiller):

  per_energy_I = {}
  per_HKL_I = {}
  for iE,Energy in enumerate(range(7090,7151)):
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
    for ikey in range(len(selected_H)):
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
  for iE,Energy in enumerate(range(7090,7151)):
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
    for ikey in range(len(selected_H)):
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

  # Now get table of intensities for the oxidized form
  # Einsle paper: Oxidized form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(III) state (absorption at higher energy, oxidized)

  # Now get table of intensities for metallic iron
  from LS49.sim.fdp_plot import george_sherrell
  class interpGS(george_sherrell):
    def fp_fdp_at_wavelength(self,angstroms):
      lookup_energy = round(12398.425/angstroms,0)
      try:
        lookup_idx = list(self.energy).index(lookup_energy)
        return self.fp[lookup_idx], self.fdp[lookup_idx]
      except Exception as e:
        lower_energies = self.energy<lookup_energy
        lower_energy = flex.max(self.energy.select(lower_energies))
        lower_idx = list(self.energy).index(lower_energy)
        higher_energies = self.energy>lookup_energy
        higher_energy = flex.min(self.energy.select(higher_energies))
        higher_idx = list(self.energy).index(higher_energy)
        #print ("breacketed by",lower_idx, lower_energy, higher_idx,higher_energy)
        # use simple linear interpolation
        energy_wt = (lookup_energy - lower_energy) / (higher_energy-lower_energy)
        fp = self.fp[lower_idx] + energy_wt * (self.fp[higher_idx] - self.fp[lower_idx])
        fdp = self.fdp[lower_idx] + energy_wt * (self.fdp[higher_idx] - self.fdp[lower_idx])
        #print (fp,fdp)
        return fp,fdp
  Fe_metallic_model = interpGS("data_sherrell/Fe.dat")

  lsq_target_function(title="Reduced form",
    label_table={"FE1":Fe_oxidized_model,"FE2":Fe_reduced_model},
    images_Gi=G.images_Gi,
    genfmodel=GF,genmiller=G)

  lsq_target_function(title="Oxidized form",
    label_table={"FE1":Fe_oxidized_model,"FE2":Fe_oxidized_model},
    images_Gi=G.images_Gi,
    genfmodel=GF,genmiller=G)

  lsq_target_function(title="Metallic iron form",
    label_table={"FE1":Fe_metallic_model,"FE2":Fe_metallic_model},
    images_Gi=G.images_Gi,
    genfmodel=GF,genmiller=G)

  exit()


  """
start here.  we are now in a position to evaluate:

pull out the slides on bandpass

DONE 1) target function for the reduced form = 2.1184329431e+12
DONE 2) target function for the oxidizd form = 2.11870886515e+12
DONE 3) target function if it is Fe-0          2.11805254818e+12
4) figure out Jacobian; figure out if I can do this by LevMar
5) with initial guess for Gi=10^-15
  """

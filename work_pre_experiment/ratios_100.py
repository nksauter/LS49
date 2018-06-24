from __future__ import print_function, absolute_import
from __future__ import division
from six.moves import range
from cctbx.array_family import flex
import pickle
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model

def generate_100():
  V = open("data.pickle","rb")
  asu = {}
  icount = 0
  while 1:
    image = pickle.load(V)
    highcc = flex.double(image["cc"]) > 0.70
    if highcc.count(True)<4: continue
    for i in range(len(image["cc"])):
      if image["cc"][i]<0.7: continue
      if asu.get(image["millers"][i],0)==0:
        icount+=1
        if icount>100: return
        yield image["millers"][i]
      asu[image["millers"][i]]=1

      #print ("CC>70%%: %20s %5.2f"%(image["millers"][i],
      #     image["cc"][i]))

if __name__=="__main__":
  M = flex.miller_index()
  for i in generate_100():
    print (i)
    M.append(i)
  W1 = 12398.425/7110.
  W2 = 12398.425/7122.

  GF = gen_fmodel(resolution=1.9,pdb_text=pdb_lines,algorithm="fft",wavelength=W1)
  GF.set_k_sol(0.435)
  #GF.make_P1_primitive()
  GF.reset_wavelength(W1)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W1)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W1)
  W1_oxidized = GF.get_intensities()
  GF.reset_wavelength(W1)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_reduced_model,newvalue=W1)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W1)
  W1_reduced = GF.get_intensities()
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_oxidized = GF.get_intensities()
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_reduced_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_reduced = GF.get_intensities()

  from cctbx import miller
  matches = miller.match_indices(W2_reduced.indices(),M)
  sel0 = flex.size_t([p[0] for p in matches.pairs()])
  sel1 = flex.size_t([p[1] for p in matches.pairs()])
  print (len(sel0))
  print (len(sel1))
  print (len(W2_reduced.indices()))

  W1_ox = W1_oxidized.select(sel0).data()
  W2_ox = W2_oxidized.select(sel0).data()
  W1_re = W1_reduced.select(sel0).data()
  W2_re = W2_reduced.select(sel0).data()
  idx   = W1_oxidized.select(sel0).indices()

  #@oxidized_ratio = W2_oxidized.select(sel0) / W1_oxidized.select(sel0)
  #reduced_ratio = W2_reduced.select(sel0) / W1_reduced.select(sel0)

  for x in range(100):
    print ("%20s %8.1f %8.1f"%(idx[x],W2_ox[x],W2_re[x]))
  ox_grid = flex.double(flex.grid(100,100))
  re_grid = flex.double(flex.grid(100,100))
  for ix in range(100):
    for iy in range(100):
      ox_grid[(ix,iy)]=W2_ox[ix]/W2_ox[iy]
      re_grid[(ix,iy)]=W2_re[ix]/W2_re[iy]
  effect = ox_grid/re_grid
  from matplotlib import pyplot as plt
  plt.imshow(effect.as_numpy_array(),cmap="bwr")
  plt.colorbar()
  plt.show()

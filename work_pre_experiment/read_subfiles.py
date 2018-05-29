from __future__ import print_function
from __future__ import division
from cctbx.array_family import flex
import pickle,glob

class GM (object):
 def generate_millers(self):
  globv = "dataX*.pickle"
  self.asu = {}
  self.icount = 0
  self.images_all = 0
  self.images_strong = 0
  for filen in glob.glob(globv):

    V = open(filen,"rb")
    while 1:
     try:
      image = pickle.load(V)
      print (image["image"])
      self.images_all+=1
      highcc = flex.double(image["cc"]) > 0.70
      if highcc.count(True)<4: continue
      self.images_strong+=1
      for i in range(len(image["cc"])):
        if image["cc"][i]<0.7: continue
        self.icount+=1
        if self.icount>10000:
            print ("readched limit")
            return
        print (filen,self.icount,"CC>70%%: %20s %5.2f"%(image["millers"][i],image["cc"][i]))
        self.asu[image["millers"][i]]=1
        yield image["millers"][i]

     except EOFError,e:
      break

if __name__=="__main__":
  M = flex.miller_index()
  G = GM()
  for i in G.generate_millers():
    M.append(i)
  print ("%d Bragg spots measured"%len(M))

  print ("%d Unique Miller indices"%(len(G.asu)))

  from LS49.sim.util_fmodel import gen_fmodel
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model

  W1 = 12398.425/7110.
  W2 = 12398.425/7122.



  GF = gen_fmodel(resolution=1.9,pdb_text=pdb_lines,algorithm="fft",wavelength=W1)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_oxidized = GF.get_intensities()
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_reduced_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  W2_reduced = GF.get_intensities()

  from cctbx import miller
  W2i = W2_reduced.indices()
  matches = miller.match_indices(M,W2i)
  sel0 = flex.size_t([p[0] for p in matches.pairs()])
  sel1 = flex.size_t([p[1] for p in matches.pairs()])
  print ("matches",len(sel0))
  sel0unique={}
  for item in sel1:
    sel0unique[W2i[item]]=1
  print ("unique",len(sel0unique))
  print ("total images %d, strong %d"%(G.images_all,G.images_strong))

  GF.xray_structure.show_summary()
  #from IPython import embed; embed()
  exit()
  W2_ox = W2_oxidized.select(sel1).data()
  W2_re = W2_reduced.select(sel1).data()
  idx   = W2_oxidized.select(sel1).indices()

  for x in xrange(100):
    print ("%20s %8.1f %8.1f"%(idx[x],W2_ox[x],W2_re[x]))
  ox_grid = flex.double(flex.grid(100,100))
  re_grid = flex.double(flex.grid(100,100))
  for ix in xrange(100):
    for iy in xrange(100):
      ox_grid[(ix,iy)]=W2_ox[ix]/W2_ox[iy]
      re_grid[(ix,iy)]=W2_re[ix]/W2_re[iy]
  effect = ox_grid/re_grid
  for x in xrange(100):
    for y in xrange(100):
      if effect[(x,y)]>1.1: effect[(x,y)]=1.1
      if effect[(x,y)]<0.9: effect[(x,y)]=0.9

  from matplotlib import pyplot as plt
  plt.imshow(effect.as_numpy_array(),cmap="bwr")
  plt.colorbar()
  plt.show()

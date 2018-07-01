from __future__ import print_function
from __future__ import division
from cctbx.array_family import flex
import glob
from six.moves import cPickle as pickle
from six.moves import range

class GM (object):
 def generate_millers(self):
  #globv = "dataX000[0-1].pickle" #test
  globv = "dataX*.pickle" # production
  self.asu = {}
  self.orig_index_Nobs = {}
  self.orig_index_Iobs = {}
  self.orig_index_millers_set = {}
  self.asu_orig_index_set = {}
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
      #highcc = flex.double(image["cc"]) > 0.80 #was 0.7 and minimum 4
      #if highcc.count(True)<3: continue
      imdict = dict()
      self.images_strong[image["image"]]=imdict
      for i in range(len(image["cc"])):
        #if image["cc"][i]<0.8: continue
        self.icount+=1
        imdict[image["millers"][i]]=dict(model=image["model"][i], obs=flex.double(image["obs"][i]),
          cc=image["cc"][i], simtbx_intensity=image["simtbx_intensity"][i],
          simtbx_miller=image["simtbx_millers"][i],orig_index=image["orig_idx"][i],
          simtbx_miller_DIALS_setting=image["simtbx_millers_DIALS_setting"][i],
        )
        #print (filen,self.icount,"CC>70%%: %20s %20s %5.2f"%(image["millers"][i],image["simtbx_millers"][i],image["cc"][i]))
        #from matplotlib import pyplot as plt
        #model = image["model"][i]
        #obs = image["obs"][i]
        #plt.plot(range(len(model)),model,"k-")
        #plt.plot(range(len(obs)),1.E10*flex.double(obs),"r-")
        #plt.show()
        self.asu[image["millers"][i]]=1


        original_index = image["orig_idx"][i]
        asu_index = image["millers"][i]
        if asu_index not in self.asu_orig_index_set:
          self.asu_orig_index_set[asu_index]=set()
        self.asu_orig_index_set[asu_index].add(original_index)

        if original_index not in self.orig_index_Nobs:
          self.orig_index_Nobs[original_index]=0
          self.orig_index_Iobs[original_index]=image["simtbx_intensity"][i]
          self.orig_index_millers_set[original_index]=set()
        else:
          assert round(self.orig_index_Iobs[original_index],5)==round(image["simtbx_intensity"][i],5)
        self.orig_index_Nobs[original_index] += 1
        self.orig_index_millers_set[original_index].add(image["millers"][i])


        yield image["millers"][i]

     except EOFError:
      break

if __name__=="__main__":
  with (open("sfall_P1_7122_amplitudes.pickle","rb")) as F:
    nominal_per_HKL_I_7122 = pickle.load(F)

  print (nominal_per_HKL_I_7122)
  nominal_per_HKL_I_7122.show_summary()

  # let's read in a P1 unit cell
  from LS49.work_for_aca_lsq.jun19_fit_Gi import get_C2_pdb_structure
  W2 = 12398.425/7122.
  GF = get_C2_pdb_structure(resolution=1.9,wavelength=W2)
  #GF.make_P1_primitive()
  uc = GF.xray_structure.unit_cell()
  print (uc)

  cb_op_C2_to_P = GF.xray_structure.change_of_basis_op_to_primitive_setting()
  nominal_C2idx = cb_op_C2_to_P.inverse().apply(nominal_per_HKL_I_7122.indices())
  nominal = {}
  for ikey in range(len(nominal_C2idx)):
    nominal[nominal_C2idx[ikey]]=(nominal_per_HKL_I_7122.data()[ikey])**2


  # confirm that all simulated data falls in the range 2.0 to 2.6 Angstrom
  # in the 2.15 to 2.45 range, how many whole-sphere Miller indices based on geometry?
  # how many do we get from original indices?
  M = flex.miller_index()
  G = GM()
  for i in G.generate_millers():
    M.append(i)
  print ("%d Bragg spots measured"%len(M))
  Norigidx = 0
  for key in G.orig_index_Nobs:  Norigidx+=G.orig_index_Nobs[key]
  print ("%d original index measurements"%Norigidx)

  print ("%d Unique Miller indices"%(len(G.asu)))
  print ("%d Unique original indices"%(len(G.orig_index_Iobs)))
  for key in sorted(G.asu_orig_index_set, key=lambda x: uc.d(x)):
    RESO = uc.d(key)
    if RESO < 2.299 or RESO > 2.301:continue
    print (uc.d(key), key,[k for k in G.asu_orig_index_set[key]],[G.orig_index_Iobs[k] for k in G.asu_orig_index_set[key]], nominal[key])
  exit("NO pickle written")
  with open("AAAAAAA.pickle","wb") as VVV:
    pickle.dump(G,VVV,pickle.HIGHEST_PROTOCOL)

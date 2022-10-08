from __future__ import division, print_function
from scitbx.array_family import flex
from libtbx.math_utils import round2
import libtbx.load_env
import os

class george_sherrell:
  def __init__(self,file):
    self.energy = flex.double()
    self.fp = flex.double()
    self.fdp = flex.double()
    with open(file,"r") as F:
      lines = F.readlines()
      for line in lines:
        tokens = [float(f) for f in line.strip().split()]
        self.energy.append(tokens[0])
        self.fp.append(tokens[1])
        self.fdp.append(tokens[2])
  def fp_fdp_at_wavelength(self,angstroms):
    lookup_energy = round2(12398.425/angstroms,0)
    lookup_idx = list(self.energy).index(lookup_energy)
    return self.fp[lookup_idx], self.fdp[lookup_idx]
  def plot_them(self,plt,f1,f2,f2a=1.0,f2b=0.0):
    plt.plot(self.energy, self.fp, f1)
    plt.plot(self.energy, f2a*self.fdp+f2b, f2)

if __name__=="__main__":

  """
  An alternate plot, using the metallic iron edge found on the web, from Exafs Materials
  of Danville, CA.  The Fe0 plot is distinctly different from the one assumed in the
  Sauter et al 2020 paper.

  Note the .dat only gives absorption; the f-prime values are 'nonsense' values and are
  not plotted here.
  Oct 8, 2022 NKS
  """
  from matplotlib import pyplot as plt
  GS = george_sherrell(os.path.join(libtbx.env.find_in_repositories("ls49_big_data"),"data_sherrell/pf-rd-ox_fftkk.out"))
  GS.plot_them(plt,f1="b*",f2="b-")
  GS = george_sherrell(os.path.join(libtbx.env.find_in_repositories("ls49_big_data"),"data_sherrell/pf-rd-red_fftkk.out"))
  GS.plot_them(plt,f1="r*",f2="r-")
  GS = george_sherrell(os.path.join(libtbx.env.find_in_repositories("ls49_big_data"),"data_sherrell/Fe.dat"))#not used here
  GS = george_sherrell(os.path.join(libtbx.env.find_in_repositories("ls49_big_data"),
                       "data_sherrell/Fe_exafs_materials_danville.dat"))
  GS.plot_them(plt,f1="m-",f2="y-",f2a=3.6,f2b=0.498)
  plt.xlim(7070, 7190)
  plt.ylim(-0.1, 5.0)

  plt.show()

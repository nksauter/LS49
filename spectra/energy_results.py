from __future__ import division, print_function
from six.moves import range
#source /reg/g/psdm/etc/psconda.sh
import numpy as np
import matplotlib.pyplot as plt

def get_results():
  import pickle
  #R = pickle.load(open("/reg/d/psdm/cxi/cxig3614/scratch/spectra209.pickle","rb"))
  R = pickle.load(open("data/spectra209.pickle","rb"))
  return R

class linear_fit:
  def __init__(self,data):
    self.x = data["expidx"]
    self.y = data["energy"]
    print(len(self.x))
    print(len(self.y))
    # y = Ap, where A = [[x 1]] and p = [[m], [c]]
    A = np.vstack([self.x, np.ones(len(self.x))]).T
    self.m,self.c = np.linalg.lstsq(A,self.y)[0]
    # y = mx + c
    # x = (1./m) y - (c/m)
  def get_residuals(self):
    print(len(self.y))
    calc_idx = (1./self.m)*np.array(self.y) - (self.c/self.m)
    print(len(calc_idx))
    return self.x - calc_idx

if __name__=="__main__":
  R = get_results()
  plt.plot(R["energy"], R["expidx"], "b.")
  plt.show()

  LF = linear_fit(R)
  e1 = 7050.
  e2 = 7150.
  idx1 = (1./LF.m)*e1 - (LF.c/LF.m)
  idx2 = (1./LF.m)*e2 - (LF.c/LF.m)

  plt.plot(R["energy"], R["expidx"], "b.")
  plt.plot([e1,e2],[idx1,idx2],"r-")
  plt.xlim([7050,7120])
  plt.show()

  residuals = LF.get_residuals()
  for image in range(len(R["energy"])):
    break
    if abs(residuals[image])>75:
      print(image,residuals[image],R["energy"][image])
      plt.plot(range(len(R['spectra'][image])),R['spectra'][image],"b-")
      plt.show()

  # plot the fitted mean energies, y = mx+c
  idx = np.array(LF.x)
  fitted_energy = LF.m*idx + LF.c
  # the histogram of the data
  n, bins, patches = plt.hist(fitted_energy, 70, range=(7050,7120),normed=1, facecolor='g', alpha=0.75)
  plt.xlabel('Energy (eV)')
  plt.ylabel('Probability')
  plt.xlim([7050,7120])
  plt.title('Histogram of Mean Pulse Energy over all Pulses')
  plt.show()

  # show a random spectrum
  spectrum_fitted_energy = LF.m * np.array(range(len(R['spectra'][200]))) + LF.c
  plt.title('Random spectrum from LG36, run 209 (event 200)')
  #plt.plot(range(len(R['spectra'][200])),R['spectra'][200],"r-")
  plt.plot(spectrum_fitted_energy,R['spectra'][200],"r-")
  plt.xlabel('Energy (eV)')
  plt.show()

  # all spectra added up
  sum_spectrum = np.array(R['spectra'][0])
  for x in range(1,100000):
    sum_spectrum += np.array(R['spectra'][x])
  plt.title('Average spectrum over LG36, run 209, 100000 events')
  #plt.plot(range(len(R['spectra'][200])),R['spectra'][200],"r-")
  plt.plot(spectrum_fitted_energy,sum_spectrum,"b-")
  plt.xlabel('Energy (eV)')
  plt.xlim([7040,7140])
  plt.show()

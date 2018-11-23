from __future__ import print_function, division
from scitbx.array_family import flex

class simple_spectrum():
  def __init__(self):
    self.eV_to_angstrom = 12398.425
    import libtbx.load_env
    import os
    self.wavelen = flex.double()
    self.response = flex.double()
    dist_dir = libtbx.env.dist_path("LS49")
    for line in open(os.path.join(dist_dir,"biocars","channel_cut_scan_wavelength.txt"),"r"):
      tokens = line.strip().split()
      if len(tokens)==0: break
      self.wavelen.append(float(tokens[0]))
      self.response.append(float(tokens[1]))
    self.energy = self.eV_to_angstrom /self.wavelen
    self.deduced_energy = flex.double(range(8000,14010,10))
    self.deduced_wavlen = self.eV_to_angstrom / self.deduced_energy
    #for x in zip(self.deduced_wavlen,self.wavelen,self.energy,self.deduced_energy,self.response):
    #  print (x)
    self.correct_baseline()
  def correct_baseline(self):
    import numpy as np
    x = self.deduced_energy[0:50].concatenate(self.deduced_energy[-50:])
    y = self.response[0:50].concatenate(self.response[-50:])
    # y = Ap, where A = [[x 1]] and p = [[m], [c]]
    A = np.vstack([x, np.ones(len(x))]).T
# workaround allows use of non-thread-safe numpy lstsq, even if openMP is enabled elsewhere in the Python program
    import os,omptbx
    workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
    omptbx.omp_set_num_threads(1)
    self.m,self.c = np.linalg.lstsq(A,y)[0]
    omptbx.omp_set_num_threads(workaround_nt)
    # y = mx + c
    # x = (1./m) y - (c/m)
    self.fit_y = float(self.m) * self.deduced_energy + float(self.c)
    self.corrected_response = self.response-self.fit_y

  def get_wavelengths(self):  return self.deduced_wavlen
  def get_flux(self):  return self.response
  def get_stats(self):
    # find median energy:
    cumulative_total = flex.sum(self.corrected_response)
    max_index = flex.max_index(self.corrected_response)
    self.mode = self.deduced_energy[max_index]
    half_cum = cumulative_total/2.
    running_cum = 0
    for idx in range(len(self.deduced_energy)):
      running_cum += self.corrected_response[idx]
      if running_cum > half_cum: break
    self.median = self.deduced_energy[idx]
    self.max_response = flex.max(self.corrected_response)
    n_channels_fwhm = (self.corrected_response > (0.5*self.max_response)).count(True)
    self.fwhm_eV = n_channels_fwhm  * 10. # 10 eV/channel
    delta_E_over_E = self.fwhm_eV / self.median
    print ("delta E / E = %.5f"%delta_E_over_E)

  def plot(self):
    from matplotlib import pyplot as plt
    #plt.plot(self.deduced_energy, self.response, "r-")
    plt.plot([self.median,self.median],[0.05,0.05+flex.max(self.corrected_response)], "b-")
    #plt.plot(self.deduced_energy, self.fit_y, "b-")
    plt.plot(self.deduced_energy, [0]*len(self.deduced_energy), "k-")
    plt.plot(self.deduced_energy, self.corrected_response,"g-")
    plt.plot([self.mode],[0.05+flex.max(self.corrected_response)], "k|")
    plt.plot([self.median-(0.5*self.fwhm_eV), self.median+(0.5*self.fwhm_eV)],[0.5*self.max_response,0.5*self.max_response], "k-")
    plt.plot([self.median-(0.5*self.fwhm_eV), self.median+(0.5*self.fwhm_eV)],[0.5*self.max_response,0.5*self.max_response], "k-")
    plt.show()

  def get_input_for_simulation(self,granularity=2):
    # decide on an energy range of 9500 eV to 13000 eV.
    selection = (self.deduced_energy >= 9500.).__and__(self.deduced_energy <= 13000.)
    staged_energy = self.deduced_energy.select(selection)
    staged_response = self.corrected_response.select(selection)
    interpolated_energy = flex.double()
    interpolated_response = flex.double()
    lenS = len(staged_energy)
    for idx in range(len(staged_energy)):
      interpolated_energy.append(staged_energy[idx])
      interpolated_response.append(staged_response[idx])
      if idx+1 == lenS: break
      for incr in range(1,int(granularity)):
        interpolated_energy.append(staged_energy[idx] + (incr/granularity)*(staged_energy[idx+1]-staged_energy[idx]))
        interpolated_response.append(staged_response[idx])

    #wavlen = self.deduced_wavlen.select(selection)
    wavlen = self.eV_to_angstrom / interpolated_energy
    #flux = self.corrected_response.select(selection)
    flux = interpolated_response
    wavelength_A = flex.min(wavlen)
    if granularity==-1:
      return wavlen[::4], flux[::4], wavelength_A
    return wavlen, flux, wavelength_A

if __name__=="__main__":
  SS = simple_spectrum()
  SS.get_stats()
  SS.plot()
  SS.get_input_for_simulation()
  print ("OK")

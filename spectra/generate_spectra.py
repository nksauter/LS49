from __future__ import division, print_function
from six.moves import range, cPickle as pickle
import numpy as np
import os

# directory location for reference files
from LS49 import ls49_big_data

def get_results():
  import six
  if six.PY3: kwargs_p = dict(encoding="latin1")
  else: kwargs_p = dict()
  R = pickle.load(open(os.path.join(ls49_big_data,"data/spectra209.pickle"),"rb"),**kwargs_p)
  return R

class linear_fit:
  def __init__(self,data):
    self.x = data["expidx"] # the expected index over the 1D spectral distribution (low pass filtered)
    self.y = data["energy"] # the ebeam energy
    #print(len(self.x))
    #print(len(self.y))
    # y = Ap, where A = [[x 1]] and p = [[m], [c]]
    A = np.vstack([self.x, np.ones(len(self.x))]).T
# workaround allows use of non-thread-safe numpy lstsq, even if openMP is enabled elsewhere in the Python program
    import os,omptbx
    workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
    omptbx.omp_set_num_threads(1)
    self.m,self.c = np.linalg.lstsq(A,self.y, rcond=-1)[0]
    # says numpy: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
    #To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
    omptbx.omp_set_num_threads(workaround_nt)
    # y = mx + c
    # x = (1./m) y - (c/m)

    # Solve the equations without numpy
    # Taken from Bevington, third edition, eqns (6.12,6.13)
    NN = len(self.x)
    from scitbx.array_family import flex
    x_array = flex.double(self.x)
    sumx = flex.sum(x_array)
    sumsq_x = sumx * sumx
    sumx_sq = flex.sum(x_array*x_array)
    sumxy = flex.sum(x_array * flex.double(self.y))
    sumy = flex.sum(flex.double(self.y))
    delta_prime = NN*sumx_sq - sumsq_x
    assert delta_prime != 0.
    mm = (1./delta_prime) * (NN*sumxy - sumx * sumy)
    cc = (1./delta_prime) * (sumx_sq*sumy - sumx * sumxy)
    from libtbx.test_utils import approx_equal
    assert approx_equal(mm,self.m)
    assert approx_equal(cc,self.c)

  def get_residuals(self):
    print(len(self.y))
    calc_idx = (1./self.m)*np.array(self.y) - (self.c/self.m)
    print(len(calc_idx))
    return self.x - calc_idx

class spectra_simulation:
  def __init__(self):
    self.R = get_results()
    self.LF = linear_fit(self.R)
    # get some information to help normalize things
    maxima = []
    bk_subtracted_sum = []
    for image in range(len(self.R["energy"])):
      bk_subtracted_sum.append(np.sum(self.R['spectra'][image]))
      maxima.append(np.max(self.R['spectra'][image]))

    self.max_of_max = max(maxima)
    average_integrated = np.mean(bk_subtracted_sum)
    print("average_integrated",average_integrated)
    self.bk_subtracted_sum = bk_subtracted_sum
    self.average_integrated = average_integrated
    self.NS = len(self.R["spectra"][0]) # number of points in each spectrum
    self.N = len(self.R["spectra"]) # number of events overall

  def plot_input_images(self,nlimit,axis="idx"):  #axis is either channel number (idx) or calibrated energy (energy)
    import matplotlib.pyplot as plt
    ylim = [-.05*self.max_of_max, 1.05*self.max_of_max]
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    for image in range(min(nlimit, self.N)):
      print(image,"ebeam = %7.2f eV"%(self.R["energy"][image]),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))
      if axis == "idx":
        plt.plot(range(self.NS),self.R['spectra'][image],"b-")
        plt.xlabel('Channel')
        plt.xlim([0,self.NS])
      elif axis == "energy":
        plt.plot(spectrum_fitted_energy,self.R['spectra'][image],"b-")
        plt.xlabel('Energy (eV)')
        ebeam, = plt.plot([self.R["energy"][image],self.R["energy"][image]],
                 [0.80*self.max_of_max, 0.90*self.max_of_max],"k-")
        expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c
        expected, = plt.plot([expected_energy,expected_energy],
                 [0.80*self.max_of_max, 0.90*self.max_of_max],"r-")
        plt.legend([ebeam,expected],["ebeam","opal"])
      plt.ylim(ylim)
      plt.title('Spectrum from LG36, run 209 (event %d)'%(image))
      plt.show()

  def plot_recast_images(self,nlimit,energy):
    import matplotlib.pyplot as plt
    ylim = [-.05*self.max_of_max, 1.05*self.max_of_max]
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset
    for image in range(min(nlimit, self.N)):
      expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
      print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))
      if True:
        plt.plot(offset_energy,self.R['spectra'][image],"b-")
        plt.xlabel('Energy (eV)')
        expected, = plt.plot([expected_energy,expected_energy],
                 [0.80*self.max_of_max, 0.90*self.max_of_max],"r-")
        plt.legend([expected],["opal"])
      plt.ylim(ylim)
      plt.title('Simulated spectrum, %.1f eV center (event %d)'%(energy,image))
      plt.show()

  def generate_recast_images(self, nlimit, energy):
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset
    for image in range(min(nlimit, self.N)):
      expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
      print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))
      yield offset_energy,self.R['spectra'][image],self.bk_subtracted_sum[image]/self.average_integrated

  def generate_recast_renormalized_images(self, nlimit, energy, total_flux):
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset
    for image in range(min(nlimit, self.N)):

      from scitbx.array_family import flex
      y = flex.double(list(self.R['spectra'][image]))
      ysum = self.bk_subtracted_sum[image]

      expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
      print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))
      assert expected_energy > 0.
      channel_flux = flex.double(100) # 100 energy channels altogether
      channel_mean_eV = flex.double(range(100)) + energy - 49.5
      eV_to_angstrom = 12398.425
      channel_wavelength = eV_to_angstrom / channel_mean_eV
      for idx in range(len(offset_energy)):
        i_energy = offset_energy[idx]
        channel = int(i_energy - (energy-50))
        if 0 <= channel < 100:
          channel_flux[channel] += self.R['spectra'][image][idx] * total_flux / self.average_integrated
      yield channel_wavelength,channel_flux,eV_to_angstrom / expected_energy

  def generate_recast_renormalized_image(self, image, energy, total_flux):
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset

    from scitbx.array_family import flex
    y = flex.double(list(self.R['spectra'][image]))
    ysum = self.bk_subtracted_sum[image]

    expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
    print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))

    channel_flux = flex.double(100) # 100 energy channels altogether
    channel_mean_eV = flex.double(range(100)) + energy - 49.5
    eV_to_angstrom = 12398.425
    channel_wavelength = eV_to_angstrom / channel_mean_eV
    for idx in range(len(offset_energy)):
        i_energy = offset_energy[idx]
        channel = int(i_energy - (energy-50))
        if 0 <= channel < 100:
          channel_flux[channel] += self.R['spectra'][image][idx] * total_flux / self.average_integrated
    yield channel_wavelength,channel_flux,eV_to_angstrom / expected_energy

  def get_average_expected_energy(self):
    idx = np.array(self.LF.x)
    fitted_energy = self.LF.m * idx + self.LF.c
    #return np.mean(fitted_energy)
    from scitbx.array_family import flex
    idx_c = flex.double(self.LF.x)
    fitted_energy_c = float(self.LF.m) * idx_c + float(self.LF.c)
    print ("numpy",np.mean(fitted_energy), "flex",flex.mean(fitted_energy_c))
    return flex.mean(fitted_energy_c)

class spectrum_simulation: # slim-down class to return only one spectrum per instance
  def __init__(self, common_data = get_results()):
    self.R = common_data # pass in get_results() please

  def select(self, image):
    self.LF = linear_fit(self.R)
    # get some information to help normalize things
    self.max_of_max = 5709559.5
    self.bk_subtracted_sum = np.sum(self.R['spectra'][image])
    self.average_integrated = 349782880.0
    self.NS = len(self.R["spectra"][0]) # number of points in each spectrum
    self.N = len(self.R["spectra"]) # number of events overall

  def plot_recast_image(self,image, energy):
    import matplotlib.pyplot as plt
    ylim = [-.05*self.max_of_max, 1.05*self.max_of_max]
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset

    expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
    print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum/self.average_integrated))
    if True:
        plt.plot(offset_energy,self.R['spectra'][image],"b-")
        plt.xlabel('Energy (eV)')
        expected, = plt.plot([expected_energy,expected_energy],
                 [0.80*self.max_of_max, 0.90*self.max_of_max],"r-")
        plt.legend([expected],["opal"])
    plt.ylim(ylim)
    plt.title('Simulated spectrum, %.1f eV center (event %d)'%(energy,image))
    plt.show()

  def generate_recast_renormalized_image(self, image, energy, total_flux):
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset

    from scitbx.array_family import flex
    y = flex.double(list(self.R['spectra'][image]))
    ysum = self.bk_subtracted_sum

    expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
    #print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
    #    self.bk_subtracted_sum/self.average_integrated))

    channel_flux = flex.double(100) # 100 energy channels altogether
    channel_mean_eV = flex.double(range(100)) + energy - 49.5
    eV_to_angstrom = 12398.425
    channel_wavelength = eV_to_angstrom / channel_mean_eV
    for idx in range(len(offset_energy)):
        i_energy = offset_energy[idx]
        channel = int(i_energy - (energy-50))
        if 0 <= channel < 100:
          channel_flux[channel] += self.R['spectra'][image][idx] * total_flux / self.average_integrated
    yield channel_wavelength,channel_flux,eV_to_angstrom / expected_energy

  def get_average_expected_energy(self):
    idx = np.array(self.LF.x)
    fitted_energy = self.LF.m * idx + self.LF.c
    from scitbx.array_family import flex
    idx_c = flex.double(self.LF.x)
    fitted_energy_c = float(self.LF.m) * idx_c + float(self.LF.c)
    return flex.mean(fitted_energy_c)

if __name__=="__main__":
  SS = spectra_simulation()
  #SS.plot_input_images(20,axis="energy")

  print("average expected energy %f eV"%(SS.get_average_expected_energy()))
  print("average ebeam is %f eB"%(np.mean(SS.R["energy"])))
  SS.plot_recast_images(20,energy=7120.)

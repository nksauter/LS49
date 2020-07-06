from __future__ import division,print_function
# %%% boilerplate specialize to packaged big data %%%
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%

from scitbx.array_family import flex
from LS49.spectra.generate_spectra import spectra_simulation
SS = spectra_simulation()
#SS.plot_recast_images(20,energy=7120.)
iterator = SS.generate_recast_renormalized_images(20,energy=7120.,total_flux=1e12)


mean_energies = flex.double()
cumulative_flux = flex.double(100)
x_axis = 12398.425/(iterator.next()[0])
N = 10000
for idx in range(N):
  iterator = SS.generate_recast_renormalized_image(image=idx,energy=7120.,total_flux=1e12)
  wavel,flux,mean = iterator.next()
  mean_energy = 12398.425 / mean
  mean_energies.append(mean_energy)
  cumulative_flux += flux
stats = flex.mean_and_variance(x_axis, cumulative_flux)
mm = stats.mean()
wsig = stats.gsl_stats_wsd()
print("Mean and standard deviation:",mm,wsig)

from matplotlib import pyplot as plt
plt.plot(x_axis, cumulative_flux/1.e12, 'r-')
plt.xlabel("energy (eV)")
plt.ylabel("photons/1 eV channel/%d images (x10^12)"%(N))
plt.title("Cumulative energy profile")
plt.axes().set_xlim((7069,7171))
plt.show()


stats = flex.mean_and_variance(mean_energies)
mu=stats.mean()
sig=stats.unweighted_sample_standard_deviation()
print (mu,sig)

plt.xlabel("energy (eV)")
plt.title("Distribution of mean pulse energies")
plt.axes().set_xlim((7069,7171))
n,bins,patches = plt.hist(mean_energies,100, range=(7070,7170), normed=0, facecolor="blue", alpha=0.75)
plt.show()

from __future__ import division, print_function
import pickle
import os
import math
os.environ["JSON_GLOB"]="null"
os.environ["PICKLE_GLOB"]="null"
os.environ["USE_POSTREFINE"]="null"
os.environ["MODEL_MODE"]="null"
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
from scitbx.array_family import flex
abc_glob_dials_refine = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_dials_refine/abcX%06d.pickle"
abc_glob_pixel_refine = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_pixel_refine/abcX%06d.pickle"
abc_glob_coarse_ground_truth = "/global/cscratch1/sd/nksauter/proj-paper1/work/abc_coverage_coarse_ground_truth/abcX%06d.pickle"
deltafast = flex.double()
deltaslow = flex.double()

from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")

from LS49.sim.util_fmodel import gen_fmodel
direct_algo_res_limit = 2.0
eV_to_angstrom = 12398.425
wavelength_A = eV_to_angstrom/7122.0

GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
GF.set_k_sol(0.435)
GF.make_P1_primitive()
sfall_main = GF.get_amplitudes()

GF.reset_wavelength(wavelength_A)
GF.reset_specific_at_wavelength(
                   label_has="FE1",tables=local_data.get("Fe_oxidized_model"),newvalue=wavelength_A)
GF.reset_specific_at_wavelength(
                   label_has="FE2",tables=local_data.get("Fe_reduced_model"),newvalue=wavelength_A)
print("USING scatterer-specific energy-dependent scattering factors")
OX = sfall_channel_oxidized = GF.get_intensities()
GF.reset_specific_at_wavelength(
                   label_has="FE1",tables=local_data.get("Fe_reduced_model"),newvalue=wavelength_A)
GF.reset_specific_at_wavelength(
                   label_has="FE2",tables=local_data.get("Fe_reduced_model"),newvalue=wavelength_A)
print("USING scatterer-specific energy-dependent scattering factors")
RE = sfall_channel_reduced = GF.get_intensities()

DD = OX.unit_cell().d(OX.indices()) # resolution of each record
SS = (DD<2.5).__and__(DD>2.1) # select on the ROI
oxss = OX.select(SS)
ress = RE.select(SS)
meanox = flex.mean(oxss.data())
diff = ress.data() - oxss.data()
qty = diff/meanox
print ("The rms difference is %10.5f%%"%(100.*math.sqrt(flex.mean(qty*qty))))
from matplotlib import pyplot as plt

fig, ax = plt.subplots()


cdf = sorted(100.*qty)
yaxis = flex.double([(y+0.5)/(len(cdf)) for y in range(len(cdf))])

from scipy.optimize import curve_fit
from scipy.special import erf
def func(xval,mu,sigma):
  return 0.5*(1.0+erf((xval-mu)/(sigma*math.sqrt(2.))))

popt, pcov = curve_fit(func, cdf, yaxis)
print ("Gauss:",popt)

ax.plot(cdf, yaxis, 'r-',label="Fig. 7 data")
ax.plot(cdf, [func(x, popt[0], popt[1]) for x in cdf], "g--", label="Gaussian")
from numpy import arctan

def Cauchyfunc(cxval,mu0,gamma):
  return 0.5 + (1./math.pi)*arctan((cxval-mu0)/gamma)
popt, pcov = curve_fit(Cauchyfunc, cdf, yaxis)
print ("Cauchy:",popt)
ax.plot(cdf, [Cauchyfunc(cx, popt[0], popt[1]) for cx in cdf], "b-",label="Lorentzian")

ax.set_xlabel('Intensity percent difference, reduced vs. oxidized (percent)')
ax.set_ylabel('CDF')
ax.set_title(r'Structure factor intensity differences')
ax.legend()
plt.show()

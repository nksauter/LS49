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
from scitbx.matrix import col
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
print ("there are %d differences"%(len(qty)))
print ("above 5%% %d"%((qty>0.05).count(True)))
print ("below -5%% %d"%((qty<-0.05).count(True)))
print ("The rms difference is %10.5f%%"%(100.*math.sqrt(flex.mean(qty*qty))))
from matplotlib import pyplot as plt

fig, ax = plt.subplots()
from matplotlib import colors
CC = plt.get_cmap("bwr")
norm = colors.Normalize(vmin=0, vmax=40)
myCC = [CC(x) for x in range(40)]

# the histogram of the data
n, bins, patches = ax.hist(100.*qty, bins=60, range=(-10.0,10.0), color= 'blue')

# fix it so the histogram is colored with a colormap
bin_centers = 0.5 * (bins[:-1] + bins[1:])
approx = (flex.double(range(60))+0.5)/60.

# scale values to interval [0,1]
col = bin_centers - min(bin_centers)
col /= max(col)

for c, p in zip(approx, patches):
    plt.setp(p, 'facecolor', CC(c))


n, bins, patches = ax.hist(100.*qty, bins=60, histtype="step",range=(-10.0,10.0), color= 'black')

ax.set_xlabel('Intensity difference, reduced vs. oxidized (percent)')
ax.set_ylabel('Count')
ax.set_title(r'Histogram of structure factor intensity differences')
plt.show()

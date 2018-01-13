from __future__ import division
from LS49.sim.util_fmodel import gen_fmodel
# generate a fake set of structure factors relevant for the step4K simulatoin
pdb_lines = open("1m2a.pdb","r").read()
GF = gen_fmodel(resolution=1.90,pdb_text=pdb_lines,algorithm="fft",wavelength=1.73)
GF.set_k_sol(0.435)
sfall_main = GF.get_amplitudes()
from scitbx.array_family import flex
amp_experimental = sfall_main.customized_copy(sigmas = flex.sqrt(sfall_main.data()))

mtz_file = "1m2a_fmodel.mtz"

mtz_out = amp_experimental.as_mtz_dataset(
      column_root_label="Fmodel",
      title="Fmodel for nanoBragg sim",
      wavelength=1.73)
mtz_obj = mtz_out.mtz_object()
mtz_obj.write(mtz_file)


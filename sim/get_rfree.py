from __future__ import division

from LS49.sim.util_fmodel import hisym_fcalc_from_pdb

pdb_lines = open("1m2a.pdb","r").read()

sfall = hisym_fcalc_from_pdb(resolution=1.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.73424)

flags = sfall.generate_r_free_flags(fraction = 0.05, max_free = 100000)

flags.show_summary(prefix = "Rfree: ")

mtz_out = sfall.as_mtz_dataset(
      column_root_label="Fcalc",
      title="Fcalc and Rfree flags",
      wavelength=1.73424)
mtz_out.add_miller_array(
      miller_array=flags,
      column_root_label="rfree_flags")
mtz_obj = mtz_out.mtz_object()
mtz_obj.write("1m2a_flags.mtz")

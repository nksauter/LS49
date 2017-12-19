from __future__ import division


def fcalc_from_pdb(resolution,pdb_text,algorithm=None,wavelength=0.9):
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None,lines = pdb_text)
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.show_summary(prefix="Input structure ")
  #
  # take a detour to insist on calculating anomalous contribution of every atom
  scatterers = xray_structure.scatterers()
  for sc in scatterers:
    from cctbx.eltbx import sasaki, henke
    #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()
  # how do we do bulk solvent?
  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  P1_primitive_xray_structure.show_summary(prefix="P1 structure ")
  fcalc = P1_primitive_xray_structure.structure_factors(
    d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
  # does provided fixed-wavelength Fcalc if anomalous_flag is True
  return fcalc.amplitudes()

def fmodel_from_pdb(resolution,pdb_text,algorithm=None,wavelength=0.9):
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None,lines = pdb_text)
  xray_structure = pdb_inp.xray_structure_simple()
  xray_structure.show_summary(prefix="Input structure ")
  #
  # take a detour to insist on calculating anomalous contribution of every atom
  scatterers = xray_structure.scatterers()
  for sc in scatterers:
    from cctbx.eltbx import sasaki, henke
    #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
    expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
    sc.fp = expected_henke.fp()
    sc.fdp = expected_henke.fdp()
  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  P1_primitive_xray_structure.show_summary(prefix="P1 structure ")
  # how do we do bulk solvent?
  import mmtbx.command_line.fmodel
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
  params2.output.type = "complex"
  params2.high_resolution = resolution
  params2.fmodel.k_sol = 0.35
  params2.fmodel.b_sol = 46.
  params2.structure_factors_accuracy.algorithm = algorithm
  import mmtbx
  f_model_complex = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = P1_primitive_xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = params2).f_model
  f_model_real = abs(f_model_complex)
  f_model_real.set_observation_type_xray_amplitude()
  f_model_real.show_summary(prefix="FMODEL ")
  return f_model_real


from __future__ import division, print_function
from scipy import constants
ENERGY_CONV = 1e10 * constants.c * constants.h / constants.electron_volt

def hisym_fcalc_from_pdb(resolution,pdb_text,algorithm=None,wavelength=0.9):
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

  fcalc = xray_structure.structure_factors(
    d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
  return fcalc.amplitudes()

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
  params2.fmodel.k_sol = 0.435 # "Babinet" optimum, used 0.35 previously
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

class gen_fmodel(object):
  def __init__(self,resolution,pdb_text,algorithm=None,wavelength=0.9):
    from iotbx import pdb
    pdb_inp = pdb.input(source_info=None,lines = pdb_text)
    xray_structure = pdb_inp.xray_structure_simple()
    xray_structure.show_summary(prefix="Input structure ")
    #
    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    from cctbx.eltbx import sasaki, henke
    for sc in scatterers:
      #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()

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

    # vvv These params restore the "legacy" solvent mask generation before
    # vvv cctbx commit 2243cc9a
    params2.mask.Fmask_res_high = 0
    params2.mask.grid_step_factor = 4
    params2.mask.solvent_radius = 1.11
    params2.mask.use_resolution_based_gridding = True
    # ^^^

    self.params2 = params2
    self.xray_structure = xray_structure
  def make_P1_primitive(self):
    primitive_xray_structure = self.xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    P1_primitive_xray_structure.show_summary(prefix="P1 structure ")
    self.cb_op_C2_to_P = self.xray_structure.change_of_basis_op_to_primitive_setting()
    self.xray_structure = P1_primitive_xray_structure
  def set_k_sol(self,newvalue):  self.params2.fmodel.k_sol = newvalue
  def reset_wavelength(self,newvalue):
    scatterers = self.xray_structure.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import sasaki, henke
      #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(newvalue)
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(newvalue)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()
  def reset_specific_at_wavelength(self,label_has,tables,newvalue,verbose=False):
    scatterers = self.xray_structure.scatterers()
    for sc in scatterers:
      if label_has in sc.label:
        newfp,newfdp = tables.fp_fdp_at_wavelength(angstroms=newvalue)
        if verbose:
          print("found",sc.element_symbol(),"label",sc.label, end=' ')
          lookup_energy=ENERGY_CONV/newvalue
          print(lookup_energy, end=' ')
          print("old",sc.fp,sc.fdp,"new",newfp,newfdp)
        sc.fp = newfp
        sc.fdp = newfdp
  def reset_specific_at_energy(self,label_has,tables,newvalue,verbose=False):
    wavelength = ENERGY_CONV/newvalue
    self.reset_specific_at_wavelength(label_has,tables,wavelength,verbose)
  def get_fmodel(self):
    import mmtbx
    return mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = self.params2)
  def get_defined_indices_fmodel(self,miller_set):
    import mmtbx
    return mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = miller_set,
      add_sigmas     = False,
      params         = self.params2)
  def get_amplitudes(self):
    import mmtbx
    f_model_complex = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = self.params2).f_model
    f_model_real = abs(f_model_complex)
    f_model_real.set_observation_type_xray_amplitude()
    return f_model_real
  def get_intensities(self):
    import mmtbx
    f_model_complex = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = self.params2).f_model
    f_model_real = f_model_complex.as_intensity_array()
    return f_model_real

from __future__ import print_function
from __future__ import division
from six.moves import cPickle as pickle
from six.moves import range
from scitbx.array_family import flex
from mmtbx.utils import fmodel_from_xray_structure
from cctbx.xray import structure_factors
from cctbx import xray
import math

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%

from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")

from LS49.sim.util_fmodel import gen_fmodel

class gen_fmodel_with_complex(gen_fmodel):
  def get_fcalc_fmodel(self):
    import mmtbx
    f_model_gen = mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = self.params2)

    self.f_calc = f_model_gen.f_obs_cache
    self.f_model = f_model_gen.f_model
    self.bulk = self.f_model.customized_copy(data = self.f_model.data() - self.f_calc.data())

  def __init__(self):
    print("Derived class")

  @classmethod
  def from_structure(cls,xray_structure,energy):
    CS = cls()
    wavelength = 12398.425/float(energy)
    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import sasaki, henke
      #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()
    CS.xray_structure = xray_structure
    CS.energy = energy
    return CS

  def from_parameters(self, resolution=1.7,
                            algorithm="fft",
                            k_sol=0.435 ):
    import mmtbx.command_line.fmodel
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    params2.output.type = "complex"
    params2.high_resolution = resolution
    params2.structure_factors_accuracy.algorithm = algorithm
    params2.fmodel.k_sol = k_sol
    params2.fmodel.b_sol = 46.
    params2.structure_factors_accuracy.grid_resolution_factor = 1/5.
    params2.mask.grid_step_factor = 10.
    self.params2 = params2
    return self

  def as_P1_primitive(self):
    primitive_xray_structure = self.xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    P1_primitive_xray_structure.show_summary(prefix="P1 structure ")
    cls = gen_fmodel_with_complex()
    cls.xray_structure = P1_primitive_xray_structure
    cls.energy = self.energy
    cls.params2 = self.params2
    cls.cb_op_C2_to_P = self.xray_structure.change_of_basis_op_to_primitive_setting()
    return cls

  def get_intensity_dictionary(self):
    W2_reduced = self.get_intensities()
    W2i = W2_reduced.indices()
    intensity_dict = {}
    for iw in range(len(W2i)):
      intensity_dict[W2_reduced.indices()[iw]] = W2_reduced.data()[iw]
    return intensity_dict


def XXXremake_intensities_at_energy(energy,FE1_model,FE2_model):

  W2 = 12398.425/float(energy)

  GF = gen_fmodel_with_complex(resolution=1.7,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.params2.fmodel.b_sol = 46.
  GF.params2.structure_factors_accuracy.grid_resolution_factor = 1/5.
  GF.params2.mask.grid_step_factor = 10.
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=FE1_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=FE2_model,newvalue=W2)
  GF.make_P1_primitive()
  GF.get_fcalc_fmodel()
  W2_reduced = GF.get_intensities()
  # Einsle paper: Reduced form has
  #    buried irons, FE1, in Fe(III) state (absorption at higher energy, oxidized)
  #    surface iron, FE2, in Fe(II) state (absorption at lower energy, reduced)

  W2i = W2_reduced.indices()
  intensity_dict = {}
  for iw in range(len(W2i)):
    intensity_dict[W2_reduced.indices()[iw]] = W2_reduced.data()[iw]
  return intensity_dict

class george_sherrell_proxy(object):
  def __init__(self,new_fp,new_fdp):
    self.fp = new_fp
    self.fdp = new_fdp
  def fp_fdp_at_wavelength(self,angstroms):
    return self.fp, self.fdp

def test_fmodel_stuff(energy,FE1_model,FE2_model):
  # section 1. C2 models.  a) whole pdb, b) no Fe, c) FE1 only, d) FE2 only
  pdb_text = local_data.get("pdb_lines")# parsing PDB structure 1M2A
  tokens = pdb_text.split("\n")
  btokens = []; ctokens = []; dtokens = []
  for token in tokens:
    splits = token.split()
    if len(splits)==12:
      if splits[11]!="FE": btokens.append(token)
    else: btokens.append(token)
    if token.find("ATOM")==0: continue
    if token.find("HETATM")==0:
      if splits[3]!="FES" or splits[2].find("FE")!=0: continue
      if splits[2]=="FE1": ctokens.append(token); continue
      if splits[2]=="FE2": dtokens.append(token); continue
    ctokens.append(token); dtokens.append(token)
  pdb_text_b = "\n".join(btokens); print (len(tokens),len(btokens),len(ctokens),len(dtokens))
  pdb_text_c = "\n".join(ctokens);
  pdb_text_d = "\n".join(dtokens);
  from iotbx import pdb
  C2_structures = []
  for tag,lines in zip(["pdb_text", "pdb_text_b", "pdb_text_c", "pdb_text_d"],
                       [pdb_text, pdb_text_b, pdb_text_c, pdb_text_d]):
    pdb_inp = pdb.input(source_info=None,lines = lines)
    xray_structure = pdb_inp.xray_structure_simple()
    xray_structure.show_summary(prefix="%s "%tag)
    C2_structures.append(xray_structure)
    print("C2 models validate OK")
  # section 2.
  if False:
    for structure in C2_structures:
     GF = gen_fmodel_with_complex.from_structure(structure,energy).from_parameters()
     ID = GF.get_intensity_dictionary()
  # section 3. Validate that the old script work2_for_aca_lsq/remake_range_intensities.py
  #            gives same intensities as the new gen_fmodel_with_complex class, for the
  #            simple case of the reduced model at energy=7125 eV
  from LS49.work2_for_aca_lsq.remake_range_intensities import \
    remake_intensities_at_energy as old_algorithm
  print("old")
  ID_ref = old_algorithm(energy,FE1_model,FE2_model)
  print("new")
  GF = gen_fmodel_with_complex.from_structure(C2_structures[0],energy).from_parameters()
  GF.reset_specific_at_energy(label_has="FE1",tables=FE1_model,newvalue=energy)
  GF.reset_specific_at_energy(label_has="FE2",tables=FE2_model,newvalue=energy)
  GF = GF.as_P1_primitive()
  ID_new = GF.get_intensity_dictionary()
  for key in ID_ref:
    assert ID_ref[key]==ID_new[key]
  assert len(ID_ref)==len(ID_new)
  print("whole structure validates OK")
  # section 4. Validate that Fmodel(whole pdb) = Fbulk + Fcalc(non-Fe) + Fcalc(FE1) + Fcalc(FE2)
  #            For this validation it is necessary from cctbx.xray import structure_factors to carry the arithmetic using A+iB, not F*F
  #GF_whole = GF #from above, but compute it again here to use algo=direct
  GF_whole = gen_fmodel_with_complex.from_structure(C2_structures[0],energy).from_parameters(algorithm="fft")
  GF_whole.reset_specific_at_energy(label_has="FE1",tables=FE1_model,newvalue=energy)
  GF_whole.reset_specific_at_energy(label_has="FE2",tables=FE2_model,newvalue=energy)
  GF_whole = GF_whole.as_P1_primitive()
  f_container = GF_whole.get_fmodel()
  Fmodel_whole = f_container.f_model
  Fcalc_whole = f_container.f_obs_cache
  Fbulk = Fcalc_whole.customized_copy(data=Fmodel_whole.data()-Fcalc_whole.data())
  GF_non_Fe = gen_fmodel_with_complex.from_structure(C2_structures[1],energy).from_parameters(algorithm="fft")
  GF_non_Fe = GF_non_Fe.as_P1_primitive()
  f_container = GF_non_Fe.get_fmodel()
  Fcalc_non_Fe = f_container.f_obs_cache
  GF_FE1 = gen_fmodel_with_complex.from_structure(C2_structures[2],energy).from_parameters(algorithm="direct")
  GF_FE1.reset_specific_at_energy(label_has="FE1",tables=FE1_model,newvalue=energy)
  GF_FE1 = GF_FE1.as_P1_primitive()
  f_container = GF_FE1.get_fmodel()
  Fcalc_FE1 = f_container.f_obs_cache
  GF_FE2 = gen_fmodel_with_complex.from_structure(C2_structures[3],energy).from_parameters(algorithm="direct")
  GF_FE2.reset_specific_at_energy(label_has="FE2",tables=FE2_model,newvalue=energy)
  GF_FE2 = GF_FE2.as_P1_primitive()
  f_container = GF_FE2.get_fmodel()
  Fcalc_FE2 = f_container.f_obs_cache
  # subpoint 4a.  Fcalc(all) = Fcalc(non-FE)+Fcalc(FE1)+Fcalc(FE2)
  test4a = Fcalc_whole.customized_copy(data = Fcalc_whole.data() -
                                              Fcalc_non_Fe.data() -
                                              Fcalc_FE1.data() - Fcalc_FE2.data())
  #print (list(test4a.data()))
  if True:
    n_outliers = 0
    for ikey,key in enumerate(test4a.indices()):
      diff0 = Fcalc_whole.data()[ikey]-Fcalc_non_Fe.data()[ikey]
      diff1 = diff0 - Fcalc_FE1.data()[ikey]
      test4a_diff = test4a.data()[ikey]
      #It's not exact because the protein Fcalc's are done by FFT algorithm, but
      # the lack of closure is generally < 1%.
      if abs(test4a_diff) > 0.01 * abs(Fcalc_whole.data()[ikey]):
        n_outliers += 1
        print ("%18s %8.2f /%8.2f %8.2f /%8.2f %8.2f /%8.2f %8.2f"%(
                               key,abs(Fcalc_whole.data()[ikey]),
                                   abs(Fcalc_non_Fe.data()[ikey]),
                                   abs(diff0),
                                   abs(Fcalc_FE1.data()[ikey]),
                                   abs(diff1),
                                   abs(Fcalc_FE2.data()[ikey]),
                                   abs(test4a_diff)
        ))
    print("Test 4a N_outliers: %d out of %d"%(n_outliers, len(test4a.indices())))
    assert (n_outliers / len(test4a.indices())) < 0.01
  # subpoint 4b.  Fmodel(whole pdb) = Fbulk + Fcalc(non-Fe) + Fcalc(FE1) + Fcalc(FE2)
  test4b = Fmodel_whole.customized_copy(data = Fmodel_whole.data() -
                                               Fbulk.data() - Fcalc_non_Fe.data() -
                                              Fcalc_FE1.data() - Fcalc_FE2.data())
  n_outliers = 0
  for ikey,key in enumerate(test4b.indices()):
      test4b_diff = test4b.data()[ikey]
      #It's not exact because the protein Fcalc's are done by FFT algorithm, but
      # the lack of closure is generally < 1%.
      if abs(test4b_diff) > 0.01 * abs(Fmodel_whole.data()[ikey]):
        n_outliers += 1
        print ("%18s %8.2f / %8.2f"%(key,abs(Fmodel_whole.data()[ikey]),
                                   abs(test4b_diff)
        ))
  print("Test 4b N_outliers: %d out of %d"%(n_outliers, len(test4b.indices())))
  assert (n_outliers / len(test4b.indices())) < 0.01

  print("Complex arithmetic validates OK")
  # section 5. Repeat the section 4 assertions, but
  # section 5a.  Assign different FE1(f',f") and FE2(f',f") & calculate Fmodel(whole)
  # section 5b.  Assemble the same answer but with explicit Python-coded Fcalcs(FE1,FE2)

  table = FE1_model
  print(table.fp_fdp_at_wavelength(angstroms = 12398.425/energy))
  new_FE1 = george_sherrell_proxy(-5,9)
  new_FE2 = george_sherrell_proxy(-4,8)
  GF_whole.reset_specific_at_energy(label_has="FE1",tables=new_FE1,newvalue=energy)
  GF_whole.reset_specific_at_energy(label_has="FE2",tables=new_FE2,newvalue=energy)
  f_container = GF_whole.get_fmodel()
  Fmodel_whole_new_fpfdp = f_container.f_model

  GF_FE1.reset_specific_at_energy(label_has="FE1",tables=new_FE1,newvalue=energy)
  f_container = GF_FE1.get_fmodel()
  Fcalc_FE1 = f_container.f_obs_cache

  GF_FE2.reset_specific_at_energy(label_has="FE2",tables=new_FE2,newvalue=energy)
  f_container = GF_FE2.get_fmodel()
  Fcalc_FE2 = f_container.f_obs_cache

  MS = Fmodel_whole_new_fpfdp.set() # avoid having to repeatedly calculate indices
  CS = Fmodel_whole_new_fpfdp.crystal_symmetry() # same here
  ALGO = structure_factors.from_scatterers(crystal_symmetry=CS,
                                           d_min=GF_whole.params2.high_resolution)
  from_scatterers_direct_fe1 = ALGO(xray_structure=GF_FE1.xray_structure,
                                    miller_set=MS,algorithm="direct")
  Fcalc_FE1_dir = from_scatterers_direct_fe1.f_calc()
  from_scatterers_direct_fe2 = ALGO(xray_structure=GF_FE2.xray_structure,
                                    miller_set=MS,algorithm="direct")
  Fcalc_FE2_dir = from_scatterers_direct_fe2.f_calc()

  for ikey,key in enumerate(MS.indices()):
        #print ("%18s %8.2f %8.2f"%(key,abs(Fcalc_FE1.data()[ikey]-Fcalc_FE1_dir.data()[ikey]),
        #                           abs(Fcalc_FE2.data()[ikey]-Fcalc_FE2_dir.data()[ikey])
        #))
    assert abs(Fcalc_FE1.data()[ikey]-Fcalc_FE1_dir.data()[ikey])==0.
    assert abs(Fcalc_FE2.data()[ikey]-Fcalc_FE2_dir.data()[ikey])==0.
  print("Test 5b, calculation with modified fp fdp (low-level interface) validates OK")

  #from IPython import embed; embed()

  """

  #xxx here we actually develop the Python code to give the A+iB for FE1 and FE2.
  looks like we should develop a new direct class to calculate F, dF/dfp, dF/dfdp
  all in one swoop, within a C++ extension module.  Use this extension to calculate
  Fcalc for FE1,FE2, instead of the conventional formalism.

  actually, gradients direct already exists. So is it possible to just wrap existing
  functions at a low-enough level to get the answer?
see tst_xray for usage.
  """

  # section 6. Validate analytical derivatives' consistency with finite differences,
  #            using correlation coefficient as a score.
  # section 6a. Get analytical derivatives from a low-level interface
  gradient_flags=xray.structure_factors.gradient_flags(
     site=False,
     u_iso=False,
     u_aniso=False,
     occupancy=False,
     fp=True,
     fdp=True)
  xray.set_scatterer_grad_flags(scatterers = GF_FE1.xray_structure.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  xray.set_scatterer_grad_flags(scatterers = GF_FE2.xray_structure.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  print ("gradients")
  sf1 = xray.ext.each_hkl_gradients_direct(
    MS.unit_cell(), MS.space_group(), MS.indices(), GF_FE1.xray_structure.scatterers(), None,
    GF_FE1.xray_structure.scattering_type_registry(), GF_FE1.xray_structure.site_symmetry_table(),
    0)
  sf2 = xray.ext.each_hkl_gradients_direct(
    MS.unit_cell(), MS.space_group(), MS.indices(), GF_FE2.xray_structure.scatterers(), None,
    GF_FE2.xray_structure.scattering_type_registry(), GF_FE2.xray_structure.site_symmetry_table(),
    0)

  #from IPython import embed; embed()
  print("Test 6a, Get analytical derivatives from a low-level interface")
  for excursion in [0.0001,0.001,0.01]:

    # finite differences
    incr_FE1_fp = george_sherrell_proxy(-5+excursion,9)
    incr_FE1_fdp = george_sherrell_proxy(-5,9+excursion)

    GF_FE1.reset_specific_at_energy(label_has="FE1",tables=incr_FE1_fp,newvalue=energy)
    f_container = GF_FE1.get_fmodel()
    incr_Fcalc_FE1_fp = f_container.f_obs_cache

    GF_FE1.reset_specific_at_energy(label_has="FE1",tables=incr_FE1_fdp,newvalue=energy)
    f_container = GF_FE1.get_fmodel()
    incr_Fcalc_FE1_fdp = f_container.f_obs_cache

    #analytical
    fp_analytical = Fcalc_FE1.data() + sf1.d_fcalc_d_fp()*excursion
    fdp_analytical = Fcalc_FE1.data() + sf1.d_fcalc_d_fdp()*excursion

    diffs_fp = fp_analytical - incr_Fcalc_FE1_fp.data()
    diffs_fdp = fdp_analytical - incr_Fcalc_FE1_fdp.data()
    #from IPython import embed; embed()

    #validation
    RMSDfp = math.sqrt( flex.sum( flex.pow( flex.abs(diffs_fp) ,2 ) ) )
    RMSDfdp = math.sqrt(flex.sum(flex.pow( flex.abs(diffs_fdp) ,2)))
    print ("with excursion %16.14f the RMSD are fp %16.14f and fdp %16.14f"%(excursion,RMSDfp,RMSDfdp))
    assert RMSDfp < 1.E-10
    assert RMSDfdp < 1.E-10


  print("Test 6b, Validate analytical vs. finite differences")

if __name__=="__main__":

  test_fmodel_stuff(energy=7125,FE1_model=Fe_oxidized_model,FE2_model=Fe_reduced_model)
  exit("OK")

  for filename,FE1,FE2 in [
    ("confirm_P1_range_reduced_intensities_dict.pickle", Fe_oxidized_model, Fe_reduced_model),
    ("confirm_P1_range_oxidized_intensities_dict.pickle", Fe_oxidized_model, Fe_oxidized_model),
    ("confirm_P1_range_metallic_intensities_dict.pickle", Fe_metallic_model, Fe_metallic_model),
    ("confirm_P1_range_swapreduced_intensities_dict.pickle", Fe_reduced_model, Fe_oxidized_model),
  ]:
    initial = remake_intensities_at_energy(7070.0,FE1,FE2)
    print ("%d keys in initial"%len(initial))

    result = {}
    for key in initial:
      result[key] = flex.double()
    print (filename)
    for incr in range(100):
      energy = 7070.5 + incr
      more = remake_intensities_at_energy(energy,FE1,FE2)
      print (energy, "with %d keys"%(len(more)))
      for key in more:
        if key in initial:
          result[key].append(more[key])
    exit("don't go here at present, dont overwrite")
    with (open(filename,"wb")) as F:
      pickle.dump(result, F, pickle.HIGHEST_PROTOCOL)

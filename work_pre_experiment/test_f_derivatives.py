from __future__ import print_function
from cctbx.array_family import flex
import pickle,glob

class gen_fmodel:
  def __init__(self,resolution,pdb_text,algorithm=None,wavelength=0.9):
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
      sc.show()

    import mmtbx.command_line.fmodel
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    params2.output.type = "complex"
    params2.high_resolution = resolution
    params2.low_resolution = 2.302
    params2.fmodel.k_sol = 0.35
    params2.fmodel.b_sol = 46.
    params2.structure_factors_accuracy.algorithm = algorithm
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
  def zero_out_specific_at_wavelength(self,label_has,verbose=False):
    scatterers = self.xray_structure.scatterers()
    for sc in scatterers:
      if label_has in sc.label:
        sc.fp = 0.000000001
        sc.fdp = 0.000000001
  def reset_specific_at_wavelength(self,label_has,tables,newvalue,verbose=False):
    scatterers = self.xray_structure.scatterers()
    for sc in scatterers:
      if label_has in sc.label:
        newfp,newfdp = tables.fp_fdp_at_wavelength(angstroms=newvalue)
        if verbose:
          lookup_energy=12398.425/newvalue
          print ("found",sc.element_symbol(),"label",sc.label,
                 lookup_energy,
                 "old",sc.fp,sc.fdp,"new",newfp,newfdp)
        sc.fp = newfp
        sc.fdp = newfdp
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
  def get_complex(self):
   import mmtbx
   f_model_complex = mmtbx.utils.fmodel_from_xray_structure(
     xray_structure = self.xray_structure,
     f_obs          = None,
     add_sigmas     = False,
     params         = self.params2)
   return f_model_complex.f_model, f_model_complex.f_calc

if __name__=="__main__":
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model

  W2 = 12398.425/7125.
  """
  # First get F's with protein, metals, and absorption
  GF = gen_fmodel(resolution=2.3,pdb_text=pdb_lines,algorithm="direct",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  ALLmodel, ALLcalc = GF.get_complex()
  I = ALLmodel.indices()
  C = ALLmodel.data()
  #print (len(I))
  #for i in xrange(len(I)):
  #  print (i,I[i],C[i])

  # Next, get F's with just protein + bulk_solvent
  plines_src = open("1m2a.pdb","r").readlines()
  plines_dst = []
  for line in plines_src:
    if "FE1" in line or "FE2" in line: continue
    plines_dst.append(line)
  pdb_lines = "".join(plines_dst)
  print (pdb_lines)
  GF = gen_fmodel(resolution=2.3,pdb_text=pdb_lines,algorithm="direct",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  Pmodel,Pcalc = GF.get_complex()

  # Next, get F's with just metals
  plines_src = open("1m2a.pdb","r").readlines()
  plines_dst = []
  for line in plines_src:
    if (not("HETATM" in line) and not("ATOM" in line)) or "FE1" in line or "FE2" in line:
      plines_dst.append(line)
  pdb_lines = "".join(plines_dst)
  print (pdb_lines)
  GF = gen_fmodel(resolution=2.3,pdb_text=pdb_lines,algorithm="direct",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  GF.reset_specific_at_wavelength(label_has="FE2",tables=Fe_reduced_model,newvalue=W2)
  Femodel,Fecalc = GF.get_complex()

  for i in xrange(len(I)):
    print ("%4d %16s%16s%16s"%(i,I[i],Fecalc.indices()[i],Pcalc.indices()[i]),Fecalc.data()[i],Pcalc.data()[i],Fecalc.data()[i]+Pcalc.data()[i],ALLcalc.data()[i])
  # shows that protein + metal = all atoms.
  # next, want to show that zeroed metal + anomalous metal = metal
  # it is not clear how to get f0(angle) function.
  # figure f0, K1, K2, ab.  Verify anomalous_metal = (f0+fp)ab
  # start working out the derivatives

  exit()
  """
  print ("Testing f0")
  # Next, get F's with just Fe1,
  plines_src = open("1m2a.pdb","r").readlines()
  plines_dst = []
  for line in plines_src:
    if (not("HETATM" in line) and not("ATOM" in line)) or ("FE1" in line and "A 201" in line):
      plines_dst.append(line)
  pdb_lines = "".join(plines_dst)
  print (pdb_lines)
  GF = gen_fmodel(resolution=1.8,pdb_text=pdb_lines,algorithm="direct",wavelength=W2)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(W2)
  GF.reset_specific_at_wavelength(label_has="FE1",tables=Fe_oxidized_model,newvalue=W2)
  Fe1model,Fe1calc = GF.get_complex()


  sc = GF.xray_structure.scatterers()[0]
  fpfdp = (complex(sc.fp,sc.fdp))

  GF.reset_wavelength(W2)
  GF.zero_out_specific_at_wavelength(label_has="FE1")
  Fe1model,Fe1zcalc = GF.get_complex()
  UC = Fe1zcalc.unit_cell()
  x = flex.double()
  y = flex.double()
  for i in xrange(len(Fe1calc.indices())):
    K1 = Fe1calc.data()[i]
    K2 = Fe1zcalc.data()[i]
    lhs = (K1/K2)-1.
    f0 = fpfdp/lhs
    f0r = complex(round(f0.real,7),round(f0.imag,7))
    f0rr = f0r.real
    print ("%4d %16s%16s"%(i,Fe1calc.indices()[i],Fe1zcalc.indices()[i]),K1,K2,lhs,f0rr)
    x.append(UC.d(Fe1calc.indices()[i]))
    y.append(f0rr)
  from matplotlib import pyplot as plt
  plt.plot(x,y,"r.")
  plt.show()

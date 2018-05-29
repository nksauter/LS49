from __future__ import print_function
from cctbx.array_family import flex
import pickle,glob

class gen_fmodel:
  def __init__(self,resolution,pdb_text,algorithm=None,wavelength=0.9,verbose=False):
    from iotbx import pdb
    pdb_inp = pdb.input(source_info=None,lines = pdb_text)
    xray_structure = pdb_inp.xray_structure_simple()
    if verbose: xray_structure.show_summary(prefix="Input structure ")
    #
    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import sasaki, henke
      #expected_sasaki = sasaki.table(sc.element_symbol()).at_angstrom(wavelength)
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()
      if verbose: sc.show()

    import mmtbx.command_line.fmodel
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    params2.output.type = "complex"
    params2.high_resolution = min(resolution)
    params2.low_resolution = max(resolution)
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
  def reset_specific_to_fpfdp(self,label_has,fp,fdp,verbose=False):
    scatterers = self.xray_structure.scatterers()
    for sc in scatterers:
      if label_has in sc.label:
        sc.fp = fp
        sc.fdp = fdp
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

def single_atom_workup(workup,resolution,angstroms):
  #print ("Testing f0 for workup %s"%(str(workup)))
  # Next, get F's with just Fe1, 
  plines_src = open("1m2a.pdb","r").readlines()
  plines_dst = []
  for line in plines_src:
    if (not("HETATM" in line) and not("ATOM" in line)) or (workup[0] in line):
      plines_dst.append(line)
  pdb_lines = "".join(plines_dst)
  #print (pdb_lines)
  GF = gen_fmodel(resolution=resolution,pdb_text=pdb_lines,algorithm="direct",
                  wavelength=angstroms,verbose=False)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(angstroms)
  GF.reset_specific_to_fpfdp(label_has=workup[1],fp=workup[2],fdp=workup[3])
  Fe1model,metalcalc = GF.get_complex()

  sc = GF.xray_structure.scatterers()[0]
  fpfdp = (complex(sc.fp,sc.fdp))   

  GF.reset_wavelength(angstroms)
  GF.zero_out_specific_at_wavelength(label_has=workup[1])
  Fe1model,metalzcalc = GF.get_complex()
  re,im = metalzcalc.data().parts()
  vec2metalzcalc = flex.vec2_double(re,im)
  vecf0 = flex.double()
  #from IPython import embed; embed()
  # workaround for lack of a mat2_double class
  #delf_delfp = flex.mat2_double()
  #delf_delfdp = flex.mat2_double()
  delf_delfp1 = flex.vec2_double()
  delf_delfdp1 = flex.vec2_double()
  delf_delfp2 = flex.vec2_double()
  delf_delfdp2 = flex.vec2_double()

  for i in xrange(len(metalcalc.indices())):
    K1 = metalcalc.data()[i]
    K2 = metalzcalc.data()[i]
    lhs = (K1/K2)-1.
    f0 = fpfdp/lhs
    f0rr = f0.real
    vecf0.append(f0rr)
    #delf_delfp.append((1./f0rr,0.,0.,1./f0rr))
    #delf_delfdp.append((0.,-1./f0rr,1./f0rr,0.))
    delf_delfp1.append((1./f0rr,0.))
    delf_delfp2.append((0.,1./f0rr))
    delf_delfdp1.append((0.,-1./f0rr,))
    delf_delfdp2.append((1./f0rr,0.))
  vec2an_bn = vec2metalzcalc
  #return delf_delfp * vec2an_bn, delf_delfdp * vec2an_bn
  return (flex.vec2_double(delf_delfp1.dot(vec2an_bn), delf_delfp2.dot(vec2an_bn)),
          flex.vec2_double(delf_delfdp1.dot(vec2an_bn), delf_delfdp2.dot(vec2an_bn)))
def util_single_atom_workup_functional_only(workup,resolution,angstroms):
  plines_src = open("1m2a.pdb","r").readlines()
  plines_dst = []
  for line in plines_src:
    if (not("HETATM" in line) and not("ATOM" in line)) or (workup[0] in line):
      plines_dst.append(line)
  pdb_lines = "".join(plines_dst)
  #print (pdb_lines)
  GF = gen_fmodel(resolution=resolution,pdb_text=pdb_lines,algorithm="direct",
                  wavelength=angstroms,verbose=False)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(angstroms)
  GF.reset_specific_to_fpfdp(label_has=workup[1],fp=workup[2],fdp=workup[3])
  Fe1model,metalcalc = GF.get_complex()
  return metalcalc
def tst_single_atom_workup():
  resolution = (2.095,2.505)
  eV = 7125.
  angstroms = eV_as_angstroms(eV)
  workup = ("FE1  FES A 201","FE1",-4.,2.0)
  FP,FDP = single_atom_workup(workup,resolution,angstroms)
  metal  = util_single_atom_workup_functional_only(workup,resolution,angstroms)
  workup1 = ("FE1  FES A 201","FE1",-3.9999,2.0)
  metal1  = util_single_atom_workup_functional_only(workup1,resolution,angstroms)
  idx = metal.indices()
  idx1 = metal1.indices()
  data = metal.data()
  data1 = metal1.data()

  from scitbx.matrix import col
  for i in range(len(idx)):
    A = (col((data1[i].real,data1[i].imag))-col((data[i].real,data[i].imag)))/0.0001
    B = col(FP[i])
    print (i,"%8.4f,%8.4f   %8.4f,%8.4f"%(A[0],A[1],B[0],B[1]))

def eV_as_angstroms(eV):
  return 12398.425/eV

def at_one_eV(eV,values,tst_delF=False):
  # values is required to be a flex double with (fp-FE1,fdp-FE1,fp-FE2,fdp-FE2)
  resolution = (2.095,2.505)
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model
  angstroms = eV_as_angstroms(eV)
  GF = gen_fmodel(resolution=resolution,pdb_text=pdb_lines,
                         algorithm="fft",wavelength=angstroms)
  GF.set_k_sol(0.435)
  GF.reset_wavelength(angstroms)
  GF.reset_specific_to_fpfdp(label_has="FE1",fp=values[0],fdp=values[1])
  GF.reset_specific_to_fpfdp(label_has="FE2",fp=values[2],fdp=values[3])
  ALLmodel, ALLcalc = GF.get_complex()
  I = ALLmodel.indices()
  All = ALLmodel.data()
  MetalloProtein = ALLcalc.data()
  Bulk = All - MetalloProtein
  re,im = All.parts()
  vec2All = flex.vec2_double(re,im)
  if tst_delF:
    delF_FE1fp = flex.vec2_double(len(I))
    delF_FE1fdp = flex.vec2_double(len(I))
    for workup in [("FE1  FES A 201","FE1",values[0],values[1]),
                   ("FE1  FES B 202","FE1",values[0],values[1])]:
      FP,FDP = single_atom_workup(workup,resolution,angstroms)
      delF_FE1fp += FP
      delF_FE1fdp += FDP
    print ("OKT")

    delF_FE2fp = flex.vec2_double(len(I))
    delF_FE2fdp = flex.vec2_double(len(I))
    for workup in [("FE2  FES A 201","FE2",values[2],values[3]),
                   ("FE2  FES B 202","FE2",values[2],values[3])]:
      FP,FDP = single_atom_workup(workup,resolution,angstroms)
      delF_FE2fp += FP
      delF_FE2fdp += FDP
    print ("OKT")
    re,im = MetalloProtein.parts()
    vec2MP = flex.vec2_double(re,im)

    return (vec2MP,delF_FE1fp, delF_FE1fdp, delF_FE2fp, delF_FE2fdp)
  # got the All data into vec2 form
  FdotdelF_FE1fp = 0.
  FdotdelF_FE1fdp = 0.
  for workup in [("FE1  FES A 201","FE1",values[0],values[1]),
                 ("FE1  FES B 202","FE1",values[0],values[1])]:
    FP,FDP = single_atom_workup(workup,resolution,angstroms)
    FdotdelF_FE1fp += vec2All.dot(FP)
    FdotdelF_FE1fdp += vec2All.dot(FDP)
    #print ("OK")

  FdotdelF_FE2fp = 0.
  FdotdelF_FE2fdp = 0.
  for workup in [("FE2  FES A 201","FE2",values[2],values[3]),
                 ("FE2  FES B 202","FE2",values[2],values[3])]:
    FP,FDP = single_atom_workup(workup,resolution,angstroms)
    FdotdelF_FE2fp += vec2All.dot(FP)
    FdotdelF_FE2fdp += vec2All.dot(FDP)
    #print ("OK")
  return (I, vec2All.dot(vec2All), FdotdelF_FE1fp, FdotdelF_FE1fdp, FdotdelF_FE2fp, FdotdelF_FE2fdp)


if __name__=="__main__":
  from LS49.sim.step5_pad import pdb_lines,Fe_oxidized_model,Fe_reduced_model
  eV = 7125.
  #values is required to be a flex double with (fp-FE1,fdp-FE1,fp-FE2,fdp-FE2)
  HASH = at_one_eV(eV=eV,values=flex.double([-4.,2.,-4.,2.]))
  print (HASH)
  # now we've taken care of F dot FP.  OK to go back to setting up an LBFGS minimizer
  # evaluating functional and gradients for 4-parameter fit.

  # the following tests the derivatives of F one FE type at a time
  HASH0 = at_one_eV(eV=eV,values=flex.double([-4.,2.,-4.,2.]),tst_delF=True)
  HASH1 = at_one_eV(eV=eV,values=flex.double([-3.9999,2.,-4.,2.]),tst_delF=True)
  HASH2 = at_one_eV(eV=eV,values=flex.double([-4.,2.0001,-4.,2.]),tst_delF=True)
  HASH3 = at_one_eV(eV=eV,values=flex.double([-4.,2.,-3.9999,2.]),tst_delF=True)
  HASH4 = at_one_eV(eV=eV,values=flex.double([-4.,2.,-4.,2.0001]),tst_delF=True)
  
  from scitbx.matrix import col
  print (HASH0)
  print (HASH1)
  # this test verifies analytical derivatives at the level of Fcalc of metalloprotein
  for i in xrange(len(HASH[0])):
    A = (col(HASH1[0][i])-col(HASH0[0][i]))/0.0001
    B = col(HASH0[1][i])
    print (i,"%8.4f,%8.4f   %8.4f,%8.4f"%(A[0],A[1],B[0],B[1]))
  # still have to test derivatives over the whole target functional
  # perform the verification with algorithm="direct", then switch to "fft" for production






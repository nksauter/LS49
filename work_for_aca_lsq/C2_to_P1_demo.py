from __future__ import division, print_function
from six.moves import range
import os, math
from mmtbx.utils import fmodel_from_xray_structure


message="""\nThe program demonstrates that f_model does not apparently adhere
to the symmetry of the model.  Structure factors computed in the
published space group (C2) should correspond exactly to those seen when
the structure is expanded to P1, except for a factor of 4 on intensities
due to the 2-fold decrease in unit cell volume. PDB structure 1M2A is
downloaded as the example.  The metric used is the r.m.s. difference in
the C2:P1 structure factor intensity ratio, away from the ideal value of 4.0,
for those Miller indices in the 2.25 - 2.35 Angstrom annulus.\n"""

class gen_fmodel:
  def __init__(self,resolution,pdb_text,wavelength=None,**kwargs):
    from iotbx import pdb
    self.d_min = resolution
    self.algorithm = kwargs["algorithm"]
    self.method=kwargs["sfall"]
    pdb_inp = pdb.input(source_info=None,lines = pdb_text)
    xray_structure = pdb_inp.xray_structure_simple()
    #xray_structure.show_summary(prefix="Input structure ")
    self.xray_structure = xray_structure

    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    for sc in scatterers:
      from cctbx.eltbx import henke
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()

    if self.method=="f_model":
      from mmtbx.programs.fmodel import master_phil as phil2
      #phil2.show()
      params2 = phil2.extract()
      # adjust the cutoff of the generated intensities to assure that
      # statistics will be reported to the desired high-resolution limit
      # even if the observed unit cell differs slightly from the reference.
      params2.high_resolution = self.d_min
      params2.output.type = "complex"
      params2.fmodel.k_sol = kwargs["k_sol"]
      params2.fmodel.b_sol = kwargs["b_sol"]
      params2.structure_factors_accuracy.algorithm = kwargs["algorithm"]
      params2.structure_factors_accuracy.grid_resolution_factor = kwargs["grid_resolution_factor"]
      if kwargs["k_sol"]>0 and kwargs["b_sol"]>0:
        params2.mask.grid_step_factor = kwargs["grid_step_factor"]
      self.params2 = params2

  def make_P1_primitive(self):
    primitive_xray_structure = self.xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    #P1_primitive_xray_structure.show_summary(prefix="P1 structure ")
    #print ("UC volume",P1_primitive_xray_structure.unit_cell().volume())
    self.cb_op_C2_to_P = self.xray_structure.change_of_basis_op_to_primitive_setting()
    self.xray_structure = P1_primitive_xray_structure

  def get_amplitudes(self):
   if self.method=="f_calc":
    f_calc_complex = self.xray_structure.structure_factors(d_min=self.d_min,
      algorithm=self.algorithm).f_calc()
    f_calc_real = abs(f_calc_complex)
    f_calc_real.set_observation_type_xray_amplitude()
    return f_calc_real

   elif self.method=="f_model":
    f_model_complex = fmodel_from_xray_structure(
      xray_structure = self.xray_structure,
      f_obs          = None,
      add_sigmas     = False,
      params         = self.params2).f_model
    f_model_real = abs(f_model_complex)
    f_model_real.set_observation_type_xray_amplitude()
    return f_model_real

def get_C2_pdb_structure(resolution,**kwargs):

    pdb_lines = open("./1m2a.pdb","r").read()
    W2 = 12398.425/7122.
    return gen_fmodel(
         resolution=resolution, pdb_text=pdb_lines,
         wavelength=W2, **kwargs
    )

def reproduce_P1_model(resolution,**kwargs):
  GF = get_C2_pdb_structure(resolution,**kwargs)
  GF.make_P1_primitive()
  sfall = GF.get_amplitudes()
  sfallf = sfall.data(); sfallidx = sfall.indices()
  #intensities = {}
  #for N in range(len(sfallidx)):
  #  intensities[sfallidx[N]]=sfallf[N] * sfallf[N]

  C2idx=GF.cb_op_C2_to_P.inverse().apply(sfallidx) # put P1 indices back into C-setting
  intensities_backtransformed_to_C2 = {}
  for N in range(len(C2idx)):
    intensities_backtransformed_to_C2[C2idx[N]]=sfallf[N] * sfallf[N]
  return intensities_backtransformed_to_C2

def reproduce_C2_model(resolution,**kwargs):
  GF = get_C2_pdb_structure(resolution,**kwargs)
  sfall = GF.get_amplitudes()
  intensities = {}
  sfallf = sfall.data(); sfallidx = sfall.indices()
  for N in range(len(sfallidx)):
    intensities[sfallidx[N]]=sfallf[N] * sfallf[N]
  return intensities

def get_unit_cell():
    from iotbx import pdb
    pdb_lines = open("./1m2a.pdb","r").read()
    pdb_inp = pdb.input(source_info=None,lines = pdb_lines)
    xray_structure = pdb_inp.xray_structure_simple()
    return xray_structure.unit_cell()

def run(**kwargs):
  # the kwargs are: sfall<f_calc|f_model>, k_sol, b_sol, algorithm<direct|fft>,
  # grid_resolution_factor <as in f_model/f_model> default 1/3.
  # grid_step_factor <as in mmtbx/masks/__init__.py> default 4.
  print()
  print(kwargs)
  C2_model = reproduce_C2_model(resolution=global_res,**kwargs)
  P1_model = reproduce_P1_model(resolution=global_res,**kwargs)

  C2_unit_cell = get_unit_cell()
  sumsq = 0.
  nrefl = 0
  for C2idx in C2_model:
    if 2.25 < C2_unit_cell.d(C2idx) < 2.35:
      nrefl += 1
      ratio = C2_model[C2idx] / (4.0*P1_model[C2idx])
      deviation = 1.-ratio
      sumsq += deviation * deviation
  print ("The rms deviation is %6.3f%%"%(100.*math.sqrt(sumsq / nrefl)))
  print()

def download_coords(struct):
  if not os.path.isfile(struct+".pdb"):
    from mmtbx.command_line.fetch_pdb import run2
    run2([struct])

if __name__=="__main__":
  print (message)
  download_coords("1m2a")
  global_res = 1.7

  run(sfall="f_calc",algorithm="direct")
  run(sfall="f_calc",algorithm="fft")

  print ("==========")
  run(sfall="f_model",k_sol=0.0,b_sol=0.0,algorithm="direct",
       grid_resolution_factor=1/3.)
  run(sfall="f_model",k_sol=0.0,b_sol=0.0,algorithm="fft",
      grid_resolution_factor=1/3.)
  run(sfall="f_model",k_sol=0.0,b_sol=0.0,algorithm="fft",
      grid_resolution_factor=1/5.)

  print ("==========")
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="direct",
      grid_resolution_factor=1/3.,grid_step_factor=4.)
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="fft",
      grid_resolution_factor=1/3.,grid_step_factor=4.)
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="fft",
      grid_resolution_factor=1/5.,grid_step_factor=4.)
  print ("==========")
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="direct",
      grid_resolution_factor=1/3.,grid_step_factor=10.)
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="fft",
      grid_resolution_factor=1/3.,grid_step_factor=10.)
  run(sfall="f_model",k_sol=0.35,b_sol=46.,algorithm="fft",
      grid_resolution_factor=1/5.,grid_step_factor=10.)

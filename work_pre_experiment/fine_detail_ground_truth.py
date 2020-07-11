from __future__ import division, absolute_import, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from cctbx import crystal_orientation
from simtbx.nanoBragg import nanoBragg
import scitbx
import math

from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49 import ls49_big_data
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step5_pad import data


def tst_all(serial_no): #emulates the action of step5_pad.py in assigning a coarse orientation to each simulation event
  Nimages = 100000
  from LS49 import legacy_random_orientations
  random_orientations = legacy_random_orientations(Nimages)
  for iteration in range(Nimages):
    rand_ori = sqr(random_orientations[iteration])
    if serial_no == iteration:
      return rand_ori
  return None

from LS49.work_pre_experiment.step5_ang_misset import get_items
class nanoBragg_mock:
 def __init__(self,mosaic_domains=25,mosaic_spread_deg=0.05):

  local_data = data()
  direct_algo_res_limit = 1.7

  wavelength_A = 1.3 # Angstroms, dummy value
  GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=local_data.get("pdb_lines"),algorithm="fft",wavelength=wavelength_A)
  GF.set_k_sol(0.435)
  GF.make_P1_primitive()
  self.sfall_main = GF.get_amplitudes()

  SIM = nanoBragg(detpixels_slowfast=(3000,3000),pixel_size_mm=0.11,Ncells_abc=(10,10,10),
    # workaround for problem with wavelength array, specify it separately in constructor.
    wavelength_A=wavelength_A,verbose=0)

  SIM.mosaic_spread_deg = mosaic_spread_deg # interpreted by UMAT_nm as a half-width stddev
  SIM.mosaic_domains = mosaic_domains       # mosaic_domains setter must come after mosaic_spread_deg setter

  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  scitbx.random.set_random_seed(1234)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=SIM.mosaic_spread_deg * math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(SIM.mosaic_domains)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )
  SIM.set_mosaic_blocks(UMAT_nm)
  self.SIM = SIM


 def set_rotation(self,rotation):

  print("Determinant",rotation.determinant())
  Amatrix_rot = (rotation * sqr(self.sfall_main.unit_cell().orthogonalization_matrix())).transpose()
  prefix = "coarse_grain"
  print("RAND_ORI", prefix)
  print(Amatrix_rot[0],Amatrix_rot[1],Amatrix_rot[2])
  print(Amatrix_rot[3],Amatrix_rot[4],Amatrix_rot[5])
  print(Amatrix_rot[6],Amatrix_rot[7],Amatrix_rot[8])
  print()

  self.SIM.Amatrix_RUB = Amatrix_rot
  #workaround for failing init_cell, use custom written Amatrix setter
  Amat = sqr(self.SIM.Amatrix).transpose() # recovered Amatrix from SIM
  from cctbx import crystal_orientation
  Ori = crystal_orientation.crystal_orientation(Amat, crystal_orientation.basis_type.reciprocal)
  print("Python unit cell from SIM state",Ori.unit_cell())
  self.MM = self.SIM.get_mosaic_domains_abc_phi_0()
 def get_average_abc(self):
   av_a = col((0.,0.,0.))
   av_b = col((0.,0.,0.))
   av_c = col((0.,0.,0.))
   for val in self.MM:
     av_a += col((val[0],val[1],val[2]))
     av_b += col((val[3],val[4],val[5]))
     av_c += col((val[6],val[7],val[8]))
   av_a/=len(self.MM); av_b/=len(self.MM); av_c/=len(self.MM);
   #further modify the average vector so it has the original length
   La = col((self.MM[0][0],self.MM[0][1],self.MM[0][2])).length()
   Lb = col((self.MM[0][3],self.MM[0][4],self.MM[0][5])).length()
   Lc = col((self.MM[0][6],self.MM[0][7],self.MM[0][8])).length()
   av_a *= La/(av_a.length())
   av_b *= Lb/(av_b.length())
   av_c *= Lc/(av_c.length())
   return (av_a[0],av_a[1],av_a[2],av_b[0],av_b[1],av_b[2],av_c[0],av_c[1],av_c[2])

if __name__ == "__main__":
  for item in get_items():
    abc = item['ABC']
    print(abc[0],abc[1],abc[2])
    print(abc[3],abc[4],abc[5])
    print(abc[6],abc[7],abc[8])
    serial_no = item['serial_no']
    rotation = tst_all(serial_no)
    M = nanoBragg_mock()
    M.set_rotation(rotation)
    print (M.get_average_abc())
    #from IPython import embed; embed()
    exit()

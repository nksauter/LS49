import os,cPickle as pickle,math
from scitbx.matrix import sqr,col
from cctbx.crystal_orientation import crystal_orientation
class ScoringContainer:
  pass

  def angular_rotation(self):
    print "model"
    model_unit = self.direct_matrix_as_unit_vectors(self.model)
    print "reference"
    reference_unit = self.direct_matrix_as_unit_vectors(self.reference,reorthogonalize=True)
    rotation = model_unit*reference_unit.inverse()
    UQ = rotation.r3_rotation_matrix_as_unit_quaternion()
    UQ = UQ.normalize() # bugfix; without this many evaluations come to 0.00000 degrees
    angle, axis = UQ.unit_quaternion_as_axis_and_angle()
    print "axis length",axis.length()
    print "axis %7.4f %7.4f %7.4f"%(axis[0],axis[1],axis[2])
    return angle

  def direct_matrix_as_unit_vectors(self,ori,reorthogonalize=False):
    direct = sqr(ori.direct_matrix())
    A = col((direct[0],direct[1],direct[2])).normalize()
    B = col((direct[3],direct[4],direct[5])).normalize()
    C = col((direct[6],direct[7],direct[8])).normalize()
    print "gamma deg",math.acos(A.dot(B))*180./math.pi
    print "alpha deg",math.acos(B.dot(C))*180./math.pi
    print " beta deg",math.acos(C.dot(A))*180./math.pi
    direct_as_unit = sqr((A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2]))
    return direct_as_unit

def image_case_factory(item_no):
  mosflm = sqr((1,0,0,0,1,0,0,0,1)) # mosflm lab vectors in their own frame
  labelit = sqr((0,0,-1,-1,0,0,0,1,0)) # mosflm basis vectors in labelit frame
  MLx=(0,0,-1)
  MLy=(-1,0,0)
  MLz=(0,1,0) # "virtual" rotation axis being used by cctbx.xfel for CXI
  LM = mosflm * labelit.inverse() # converts labelit frame coords to mosflm frame
  SWAPXY = sqr((0,1,0,1,0,0,0,0,1)) # in labelit frame
  SWAPZ  = sqr((1,0,0,0,1,0,0,0,-1)) # in labelit frame
  R90    = sqr((0,-1,0,1,0,0,0,0,1)) # in labelit frame, rotation 90 on beam axis

  CONTAINER_SZ=1000
  WAVELENGTH = 1.32 # given by James Holton

  container_no = item//CONTAINER_SZ
  filename = "/net/viper/raid1/sauter/fake/result/basic00/"
  filename = os.path.join( filename, "r%02d"%container_no, TRIAL,
                           "integration", "int-data_%05d.pickle"%item_no)
  trial_results = pickle.load(open(filename,"rb"))

  current_hexagonal_ori = trial_results["current_orientation"][0]
  current_cb_op_to_primitive = trial_results["current_cb_op_to_primitive"][0]
  current_triclinic_ori = current_hexagonal_ori.change_basis(current_cb_op_to_primitive)
  filename = "/net/viper/raid1/sauter/fake/holton/mosflm_matrix"
  filename = os.path.join( filename, "%02d"%container_no, "%05d.mat"%item_no)
  lines = open(filename).readlines()

  A0 = lines[0].strip().split()
  A1 = lines[1].strip().split()
  A2 = lines[2].strip().split()
  A = sqr((float(A0[0]), float(A0[1]), float(A0[2]),
           float(A1[0]), float(A1[1]), float(A1[2]),
           float(A2[0]), float(A2[1]), float(A2[2])))

  A = A/WAVELENGTH
  Holton_hexagonal_ori = crystal_orientation(SWAPZ*SWAPXY*R90*LM*A,True)
  Holton_triclinic_ori = Holton_hexagonal_ori.change_basis(current_cb_op_to_primitive)

  c_inv_r_best = Holton_triclinic_ori.best_similarity_transformation(
               other=current_triclinic_ori,
               fractional_length_tolerance=50.,
               unimodular_generator_range=1)
  c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
  from cctbx import sgtbx
  c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
  cb_op = sgtbx.change_of_basis_op(c_inv)
  comparison_triclinic = Holton_triclinic_ori.change_basis(cb_op)
  comparison_hexagonal = comparison_triclinic.change_basis(current_cb_op_to_primitive.inverse())
  print sqr(current_hexagonal_ori.direct_matrix())-sqr(comparison_hexagonal.direct_matrix())
  print "item %d"%item_no

  SC = ScoringContainer()
  SC.model = current_hexagonal_ori
  SC.reference = comparison_hexagonal
  return SC

if __name__=="__main__":
  import sys
  TRIAL = sys.argv[1] # like "006"
  from scitbx.array_family import flex
  angular = flex.double()
  for item in xrange(1,20001):
    try:
      SC = image_case_factory(item)
      angle_deg = SC.angular_rotation()*180./math.pi
      angular.append(angle_deg)
      print "Item %5d angular offset is %8.5f deg."%(item,angle_deg)
    except IOError:
      pass
    except RuntimeError:
      print "item %d couldn't find good orientational fit."%item
    #break
  print "RMSD", math.sqrt(flex.mean(angular*angular))

from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
import libtbx.load_env # possibly implicit
import math
from libtbx.test_utils import approx_equal

class plotter:
  def __init__(self, tophat, normal):
    # take U-mats from two different distributions, apply them to unit vectors, and plot
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(2, 3,figsize=(24,14))

    # columns plot the transformation of x, y, and z unit vectors
    rows = [tophat,normal]
    for irow,dist in enumerate(rows):
      iaxes = axes[irow];
      for icol, permutation in enumerate([(0,1,2), (2,0,1), (1,2,0)]):
        axis = iaxes[icol]
        unit = self.permute_vector(vector=(1,0,0), perm = permutation)
        perm_y = self.permute_vector(vector=(0,1,0), perm = permutation).index(1)
        perm_z = self.permute_vector(vector=(0,0,1), perm = permutation).index(1)
        a2 = flex.double(); a3 = flex.double()
        for u in dist:
          U = sqr(u)
          newvec = U * unit
          #print newvec.elems
          a2.append( (180./math.pi) * math.asin(newvec[perm_y]) )
          a3.append( (180./math.pi) * math.asin(newvec[perm_z]) ) # now in degrees
        form = {0:'r.', 1:'b.'}[irow]
        axis.plot (a2,a3,form)
        axis.set_aspect("equal")
        axis.set_title("Transformation of unit vector %s"%(str(unit)))
        axis.set_xlim(-0.11,0.11)
        axis.set_ylim(-0.11,0.11)

    plt.show()

  def permute_vector(self, vector, perm):
    return (vector[perm[0]], vector[perm[1]], vector[perm[2]],)

class plotter2:  # compare the transformation of 001 with that of .57735,.57735,.57735
  def __init__(self, tophat, normal, plot):
    # take U-mats from two different distributions, apply them to unit vectors, and plot
    if plot:
      from matplotlib import pyplot as plt
      fig, axes = plt.subplots(2, 2,figsize=(8,7))
    else:
      axes = ((1,2),(3,4)) #dummy

    # columns plot the transformation of x, y, and z unit vectors
    rows = [tophat,normal]
    differences = []
    for irow,dist in enumerate(rows):
      iaxes = axes[irow];
      cube_diag = math.sqrt(1./3) # 0.57735
      for icol, RLP in enumerate([(0,0,1), (cube_diag, cube_diag, cube_diag)]):
        RLP = col(RLP)
        print "(%7.5f %7.5f %7.5f)"%(RLP.elems),"Vector length:%8.6f"%(RLP.length()),
        axis = iaxes[icol]
        unit = RLP.normalize()
        seed = col((1,0,0))
        perm2 = unit.cross(seed)
        perm3 = unit.cross(perm2)
        a2 = flex.double(); a3 = flex.double()
        difference_vectors = flex.vec3_double()
        for u in dist:
          U = sqr(u)
          newvec = U * RLP
          difference_vectors.append( newvec-RLP )
          a2.append( (180./math.pi) * math.asin(newvec.dot(perm2)) )
          a3.append( (180./math.pi) * math.asin(newvec.dot(perm3)) ) # now in degrees
        rms = math.sqrt( flex.mean ( difference_vectors.dot(difference_vectors) ) )
        print "The rms difference is", rms
        differences.append(rms)
        if plot:
          axis.plot (a2,a3,'r.')
          axis.set_aspect("equal")
          axis.set_title("Transformation of vector %s"%(str(RLP.elems)))
          axis.set_xlim(-0.11,0.11)
          axis.set_ylim(-0.11,0.11)

      assert approx_equal(differences[0],differences[1],eps=1e-04), \
      "RMS mosaic distribution for axis vector and diagonal vector should be similar, as proposed by J Holton"

    if plot: plt.show()

class check_distributions:
  @staticmethod
  def get_angular_rotation(dist):
    angle_deg = flex.double()
    for umat in dist:
      angle,axis = sqr(
        umat).r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
      angle_deg.append(angle)
    return angle_deg

MOSAIC_SPREAD = 0.05 # top hat half width rotation in degrees
SAMPLE_SIZE = 90000

def run_sim2smv(fileout):
    serial_no = 1 # look at image serial number
    from LS49.work_pre_experiment.fine_detail_ground_truth import tst_all, nanoBragg_mock
    rotation = tst_all(serial_no)
    M = nanoBragg_mock(mosaic_domains=25)
    M.set_rotation(rotation)
    UMAT_th = M.SIM.get_mosaic_blocks() # extract actually-used U-mats

    #sanity checks on the top hat distribution
    th_angles = check_distributions.get_angular_rotation(UMAT_th)
    max_angle = flex.max(th_angles)
    # not True for a Gaussian! assert max_angle <= MOSAIC_SPREAD + 0.0000001 # need to allow a small epsilon
    assert max_angle > 0.99 * MOSAIC_SPREAD # insist that max angle is near the limit we gave it
    rms_angle = math.sqrt(flex.mean(th_angles*th_angles))
    print rms_angle
    # assert rms_angle < MOSAIC_SPREAD

    M = nanoBragg_mock(mosaic_domains=200)
    M.set_rotation(rotation)
    UMAT_nm = M.SIM.get_mosaic_blocks() # extract actually-used U-mats

    #sanity check on the gaussian distribution
    nm_angles = check_distributions.get_angular_rotation(UMAT_nm)
    nm_rms_angle = math.sqrt(flex.mean(nm_angles*nm_angles))
    print nm_rms_angle
    assert approx_equal(rms_angle,nm_rms_angle,eps=1e-01), \
    "The top hat and gaussian models should have similar standard deviations"

    return UMAT_th, UMAT_nm

def tst_all(make_plots):
  fileout = "dummy_file.cbf"
  UM_th, UM_nm = run_sim2smv(fileout)
  if False and make_plots:
    P = plotter(UM_th,UM_nm)
  Q = plotter2(UM_th,UM_nm,make_plots) # suggested by Holton:
  """apply all the UMATs to a particular starting unit vector and take the rms of the resulting end points.
    You should not see a difference between starting with a vector with x,y,z = 0,0,1 vs
    x,y,z = 0.57735,0.57735,0.57735.  But if the implementation is wrong and the UMATs are being made by
    generating Gaussian-random rotations about the three principle axes, then the variance of
    0,0,1 will be significantly smaller than that of 0.57735,0.57735,0.57735."""

if __name__=="__main__":
  import sys
  make_plots = "--plot" in sys.argv
  tst_all(make_plots)
  print "OK"

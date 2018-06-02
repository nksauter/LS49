from __future__ import division
from scitbx.array_family import flex # implicit import
from scitbx.matrix import sqr,col
from xfel.command_line.print_pickle import generate_data_from_streams
from dxtbx.format.Registry import Registry
from cctbx import crystal_orientation
from dxtbx.model import Crystal
import glob,math

# some parameters
pickle_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/int*.pickle"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
global format_class

def get_coord_shifts():
  for a1 in range(-1,2):
   for a2 in range(-1,2):
    for a3 in range(-1,2):
     for a4 in range(-1,2):
      for a5 in range(-1,2):
       for a6 in range(-1,2):
        for a7 in range(-1,2):
          for a8 in range(-1,2):
           for a9 in range(-1,2):
             elems = [a1,a2,a3,a4,a5,a6,a7,a8,a9]
             mat3 = sqr(elems)
             if mat3.determinant() in [-1,1]:
               yield mat3

def get_items():
  file_list = glob.glob(pickle_glob)
  format_class = None
  for item in file_list:
    #print item
    serial_no = int(item[-16:-11])
    #print serial_no
    image_file = image_glob%serial_no
    #print image_file
    if format_class is None:
      format_class = Registry.find(image_file)
    i = format_class(image_file)
    Z = i.get_smv_header(image_file)
    ABC = Z[1]["DIRECT_SPACE_ABC"]
    abc = tuple([float(a) for a in ABC.split(",")])
    G = generate_data_from_streams([item])
    stills_process = next(G)
    yield(dict(serial_no=serial_no,ABC=abc,stills_process=stills_process))
    #for key in stills_process:
    #  print key, stills_process[key]
    #exit()

def plot_unit_cell(ax,Ori):
  #from IPython import embed; embed()
  def pvec(P,Q,c):
    ax.plot([P[0],Q[0]],[P[1],Q[1]],[P[2],Q[2]],c=c)

  direct = Ori.direct_matrix()
  zero=col((0,0,0))
  veca = col((direct[0],direct[1],direct[2]))
  vecb = col((direct[3],direct[4],direct[5]))
  vecc = col((direct[6],direct[7],direct[8]))
  print "Lengths %9.3f %9.3f %9.3f"%(veca.length(), vecb.length(),vecc.length())
  pvec(zero,veca,'r')
  pvec(zero,vecb,'g')
  pvec(zero,vecc,'b')
  pvec(veca,veca+vecb,'g')
  pvec(veca,veca+vecc,'b')
  pvec(vecb,vecb+veca,'r')
  pvec(vecb,vecb+vecc,'b')
  pvec(vecc,vecc+veca,'r')
  pvec(vecc,vecc+vecb,'g')
  pvec(veca+vecb,veca+vecb+vecc,'b')
  pvec(veca+vecc,veca+vecb+vecc,'g')
  pvec(vecc+vecb,veca+vecb+vecc,'r')


if __name__=="__main__":
  #for iu,uni in enumerate(get_coord_shifts()):
  #  print iu, uni.determinant()

  pdb_lines = open("/net/dials/raid1/sauter/LS49/1m2a.pdb","r").read()
  from LS49.sim.util_fmodel import gen_fmodel
  GF = gen_fmodel(resolution=3.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
  CB_OP_C = GF.xray_structure.change_of_basis_op_to_primitive_setting() # from C to P
  mosflm = sqr((1,0,0,0,1,0,0,0,1)) # mosflm lab vectors in their own frame
  labelit = sqr((0,0,-1,-1,0,0,0,1,0)) # mosflm basis vectors in labelit frame
  MLx=(0,0,-1)
  MLy=(-1,0,0)
  MLz=(0,1,0) # "virtual" rotation axis being used by cctbx.xfel for CXI
  LM = mosflm * labelit.inverse() # converts labelit frame coords to mosflm frame
  SWAPXY = sqr((0,1,0,1,0,0,0,0,1)) # in labelit frame
  SWAPZ  = sqr((1,0,0,0,1,0,0,0,-1)) # in labelit frame
  R90    = sqr((0,-1,0,1,0,0,0,0,1)) # in labelit frame, rotation 90 on beam axis
  LI = sqr((0,0,1,1,0,0,0,1,0))

  for stuff in get_items():
   for testuni in get_coord_shifts():
    print stuff
    print
    print stuff["ABC"]
    dx_cryst = Crystal(real_space_a=stuff["ABC"][0:3],
                       real_space_b=stuff["ABC"][3:6],
                       real_space_c=stuff["ABC"][6:9],
                       space_group_symbol="C 2")
    Ori = crystal_orientation.crystal_orientation(
      testuni*sqr(stuff["ABC"]), crystal_orientation.basis_type.direct)
    print Ori.unit_cell().parameters()
    #from IPython import embed; embed()
    sim_U = Ori.crystal_rotation_matrix()

    refined_ori = stuff["stills_process"]["current_orientation"][0]
    current_cb_op_to_primitive = stuff["stills_process"]["current_cb_op_to_primitive"][0]
    refined_triclinic_ori = refined_ori.change_basis(current_cb_op_to_primitive)
    print "UNI",testuni,sqr(testuni).determinant()
    ground_truth_ori = Ori
    print "refined ori--> %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"%(
    refined_triclinic_ori.unit_cell().parameters())
    print "truth ori  --> %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"%(
    ground_truth_ori.unit_cell().parameters())

    gamma = ground_truth_ori.unit_cell().parameters()[5]
    alpha = ground_truth_ori.unit_cell().parameters()[3]
    beta  = ground_truth_ori.unit_cell().parameters()[4]
    if gamma<80 or gamma>85 or alpha<73 or alpha>77 or beta<73 or beta>77: continue

    c_inv_r_best = ground_truth_ori.best_similarity_transformation(
               other=refined_triclinic_ori,
               fractional_length_tolerance=50.,
               unimodular_generator_range=1)
    print (sqr(c_inv_r_best)).determinant()
    #c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
    c_inv_r_int = tuple([int(round(ij,0)) for ij in c_inv_r_best])
    from cctbx import sgtbx
    c_inv = sgtbx.rt_mx(sgtbx.rot_mx(c_inv_r_int))
    cb_op = sgtbx.change_of_basis_op(c_inv)
    print cb_op;
    comparison_triclinic = ground_truth_ori.change_basis(cb_op)
    print "comparison---> %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"%(
    comparison_triclinic.unit_cell().parameters())

    missetting_rot = sqr(refined_ori.direct_matrix()) * sqr(ground_truth_ori.direct_matrix()).inverse()
    print "deter",missetting_rot.determinant()
    from LS49.missetting.mark0 import ScoringContainer
    SC = ScoringContainer()
    SC.model = refined_triclinic_ori
    SC.reference = ground_truth_ori
    angle_deg = SC.angular_rotation()*180./math.pi
    #print missetting_rot.elems
    angle,axis = missetting_rot.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
    print "Angle %8.5f degrees"%angle_deg
    print angle
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D # import dependency
    import numpy as np
    import matplotlib.pyplot as plt
    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure(figsize=(10, 8))

    #fig.suptitle('{} randomly selected cells out of total {} images'
    #             ''.format(1,2), fontsize=18)

    #ax = fig.add_subplot(111, projection='3d')
    ax = fig.gca( projection='3d')
    unit = np.linspace(0,10,100)
    zeros = 0.0 * unit
    ax.plot(unit,zeros,zeros,c='r',label="x-axis")
    ax.legend()
    ax.plot(zeros,unit,zeros,c='g',label="y-axis")
    ax.legend()
    ax.plot(zeros,zeros,unit,c='b',label="z-axis")
    ax.legend()
    test = R90*LI*LI
    print "Test",test.determinant()
    trial_ori = ground_truth_ori.change_basis(test)

    #plot_unit_cell(ax,comparison_triclinic)
    plot_unit_cell(ax,refined_triclinic_ori)
    plot_unit_cell(ax,ground_truth_ori)

    ax.set_aspect("equal")
    plt.show()
    #from IPython import embed; embed()
    #exit()
   break
# try to think about this in terms of round-tripping.
# cubic...orthorhombic...monoclinic...triclinic
# get single script to do all, thereby testing image_simulator and indexing all together
# start with a Crystal (A matrix), simulate image, index and verify A matrix

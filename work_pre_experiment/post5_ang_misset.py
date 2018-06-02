from __future__ import division, absolute_import, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr
from dxtbx.format.Registry import Registry
from cctbx import crystal_orientation

import glob,math
import dials

# some parameters
json_glob = "/net/dials/raid1/sauter/LS49_integ_step5cori/idx*.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
global format_class

def get_items():
  file_list = glob.glob(json_glob)
  format_class = None
  for item in file_list:
    serial_no = int(item[-37:-32])
    image_file = image_glob%serial_no
    #print image_file
    if format_class is None:
      format_class = Registry.find(image_file)
    i = format_class(image_file)
    Z = i.get_smv_header(image_file)
    ABC = Z[1]["DIRECT_SPACE_ABC"]
    abc = tuple([float(a) for a in ABC.split(",")])

    from dxtbx.model.experiment_list import ExperimentListFactory
    EC = ExperimentListFactory.from_json_file(item,check_format=False)[0].crystal
    try:
      yield(dict(serial_no=serial_no,ABC=abc,integrated_crystal_model=EC,postref=postrefined[serial_no]))
    except KeyError:
      continue

def parse_postrefine():
  lines = open("/net/dials/raid1/sauter/LS49_merge/merge5_redo2.log")
  result = {}
  for line in lines:
    if "ASTAR" not in line: continue
    stripped = line.replace("(","").replace(")","").replace(",","")
    tokens = stripped.split()
    try:
      serial_no = int(tokens[1][-17:-11])
      rotmat = sqr(tuple([float(a) for a in tokens[2:]]))
      assert len(rotmat)==9
      result[serial_no] = rotmat
    except (IndexError,ValueError,AssertionError):
      continue
  return result
  # a lot of files have dropped out because of thread-to-thread contention in the log
  # file (thousands).   Should re-run this sometime with nproc=1

if __name__=="__main__":
  postrefined = parse_postrefine()
  print(len(postrefined),"postrefined files")

  pdb_lines = open("/net/dials/raid1/sauter/LS49/1m2a.pdb","r").read()
  from LS49.sim.util_fmodel import gen_fmodel
  GF = gen_fmodel(resolution=3.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.7)
  CB_OP_C_P = GF.xray_structure.change_of_basis_op_to_primitive_setting() # from C to P
  print(str(CB_OP_C_P))

  icount=0
  from scitbx.array_family import flex
  angles=flex.double()
  for stuff in get_items():
    #print stuff
    icount+=1
    print("Iteration",icount)
    # work up the crystal model from integration
    direct_A = stuff["integrated_crystal_model"].get_A_inverse_as_sqr()
    permute = sqr((0,0,1,0,1,0,-1,0,0))
    sim_compatible = direct_A*permute # permute columns when post multiplying
    from cctbx import crystal_orientation
    integrated_Ori = crystal_orientation.crystal_orientation(sim_compatible, crystal_orientation.basis_type.direct)
    #integrated_Ori.show(legend="integrated")

    # work up the crystal model from postrefinement
    direct_A = stuff["postref"].inverse()
    permute = sqr((0,0,1,0,1,0,-1,0,0))
    sim_compatible = direct_A*permute # permute columns when post multiplying
    from cctbx import crystal_orientation
    postref_Ori = crystal_orientation.crystal_orientation(sim_compatible, crystal_orientation.basis_type.direct)

    # work up the ground truth from header
    header_Ori = crystal_orientation.crystal_orientation(stuff["ABC"], crystal_orientation.basis_type.direct)
    #header_Ori.show(legend="header_Ori")

    C2_ground_truth = header_Ori.change_basis(CB_OP_C_P.inverse())
    C2_ground_truth.show(legend="C2_ground_truth")

    # align integrated model with ground truth
    cb_op_align = integrated_Ori.best_similarity_transformation(C2_ground_truth,50,1)
    aligned_Ori = integrated_Ori.change_basis(sqr(cb_op_align))
    aligned_Ori.show(legend="integrated, aligned")
    print("alignment matrix", cb_op_align)

    # align postref model with ground truth
    cb_op_align = postref_Ori.best_similarity_transformation(C2_ground_truth,50,1)
    alipost_Ori = postref_Ori.change_basis(sqr(cb_op_align))
    alipost_Ori.show(legend="postrefined, aligned")
    print("alignment matrix", cb_op_align)

    from libtbx.test_utils import approx_equal
    U_postref = alipost_Ori.get_U_as_sqr()
    assert approx_equal(U_postref.determinant(),1.0)
    U_integrated = aligned_Ori.get_U_as_sqr()
    assert approx_equal(U_integrated.determinant(),1.0)
    U_ground_truth = C2_ground_truth.get_U_as_sqr()
    assert approx_equal(U_ground_truth.determinant(),1.0)

    missetting_rot = U_postref * U_ground_truth.inverse()

    angle,axis = missetting_rot.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)

    # now calculate the angle as mean a_to_a,b_to_b,c_to_c
    aoff = alipost_Ori.a.angle(C2_ground_truth.a,deg=True)
    boff = alipost_Ori.b.angle(C2_ground_truth.b,deg=True)
    coff = alipost_Ori.c.angle(C2_ground_truth.c,deg=True)
    # solved:  the reason missetting_rot doesn't exactly align postref and ground_truth is
    # that it's a monoclinic lattice, not orthorhombic.  Difference in the beta angle prevents exact alignment

    hyp = flex.mean(flex.double((aoff,boff,coff)))

    angles.append(hyp)
    print("Item %5d angular offset is %12.9f deg."%(icount,hyp))
  print("RMSD", math.sqrt(flex.mean(angles*angles)))

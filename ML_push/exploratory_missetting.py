from __future__ import division, print_function
import math
from scitbx.matrix import sqr
from LS49.work_pre_experiment.post5_ang_misset import get_items
from LS49.work_pre_experiment.post5_ang_misset import parse_postrefine # from postrefinement log "ASTAR"
# need to give the postrefinement log (demangled) as sys.argv[1]

# just exploratory.  Done in the C2 setting, rather than the P1 setting used by abc_coverage.  Return to that later.
def metric(testOri,refOri):
    # now calculate the angle as mean a_to_a,b_to_b,c_to_c
    aoff = testOri.a.angle(refOri.a,deg=True)
    boff = testOri.b.angle(refOri.b,deg=True)
    coff = testOri.c.angle(refOri.c,deg=True)
    # solved:  the reason missetting_rot doesn't exactly align postref and ground_truth is
    # that it's a monoclinic lattice, not orthorhombic.  Difference in the beta angle prevents exact alignment
    hyp = flex.mean(flex.double((aoff,boff,coff)))
    return hyp

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
  for stuff in get_items(postrefined):
    #print(stuff)
    icount+=1
    if icount>10000: break
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
    #print("alignment matrix", cb_op_align)
    metric_val = metric(aligned_Ori,C2_ground_truth)
    print("Item %5d image_no %d drefinery angular offset is %12.9f deg."%(icount,stuff["serial_no"], metric_val))

    minimum = metric_val
    #for ix in xrange(-10,11):
    #  xrotated_Ori = integrated_Ori.rotate_thru((0,0,1),ix*0.01*math.pi/180.)
    #  for iy in xrange(-10,11):
    #    yrotated_Ori = xrotated_Ori.rotate_thru((0,1,0),iy*0.01*math.pi/180.)
    #    new_aligned_ori = yrotated_Ori.change_basis(sqr(cb_op_align))
    #    grid_metric_val = metric(new_aligned_ori,C2_ground_truth)
    #    print("ix %4d"%ix,"iy %4d"%iy,"grid search angular offset is %12.9f deg."%(grid_metric_val))
    #    minimum = min(minimum, grid_metric_val)
    for ix in xrange(-10,11):
      xrotated_Ori = integrated_Ori.rotate_thru((0,0,1),ix*0.01*math.pi/180.)
      for iy in xrange(-10,11):
        yrotated_Ori = xrotated_Ori.rotate_thru((0,1,0),iy*0.01*math.pi/180.)
        for iz in xrange(-10,11):
          zrotated_Ori = yrotated_Ori.rotate_thru((1,0,0),iz*0.01*math.pi/180.)
          new_aligned_ori = zrotated_Ori.change_basis(sqr(cb_op_align))
          grid_metric_val = metric(new_aligned_ori,C2_ground_truth)
          minimum = min(minimum, grid_metric_val)
    print("the minimum angular offset is %12.9f deg"%minimum)

    # align postref model with ground truth
    cb_op_align = postref_Ori.best_similarity_transformation(C2_ground_truth,50,1)
    alipost_Ori = postref_Ori.change_basis(sqr(cb_op_align))
    alipost_Ori.show(legend="postrefined, aligned")
    #print("alignment matrix", cb_op_align)
    print("Item %5d image_no %d rs_hybrid angular offset is %12.9f deg."%(icount,stuff["serial_no"], metric(alipost_Ori,C2_ground_truth)))

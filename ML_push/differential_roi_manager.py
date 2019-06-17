from __future__ import print_function, division
from six.moves import range
from scitbx.array_family import flex
from LS49.sim.util_fmodel import gen_fmodel
from LS49.sim.step5_pad import data
from cctbx import crystal_orientation
from scitbx.matrix import sqr

def ersatz_all_orientations(N_total=100000):
    ori_N_total = N_total # number of items to simulate
    mt = flex.mersenne_twister(seed=0)
    random_orientations = []
    for iteration in range(ori_N_total):
      random_orientations.append( mt.random_double_r3_rotation_matrix() )
    return random_orientations

class differential_roi_manager(object):
  def __init__(self,idx):
    self.gen_fmodel_adapt()
    self.models4 = self.get_idx_rotation_models(idx)
    for model in self.models4:
      for roi in group of spots:
        perform util_partiality run_sim2smv
        display results

  def gen_fmodel_adapt(self):
    direct_algo_res_limit = 1.7
    self.GF = gen_fmodel(resolution=direct_algo_res_limit,pdb_text=data(
         ).get("pdb_lines"),algorithm="fft",wavelength=7122)
    self.CB_OP_C_P = self.GF.xray_structure.change_of_basis_op_to_primitive_setting() # from C to P, for ersatz model
    self.GF.set_k_sol(0.435)
    self.GF.make_P1_primitive()
    self.sfall_main = self.GF.get_amplitudes()

  def get_idx_rotation_models(self,idx):
    rotation = sqr(ersatz_all_orientations()[idx])
    Amat = (rotation * sqr(self.sfall_main.unit_cell().orthogonalization_matrix())).transpose()
    from LS49.work_pre_experiment.post5_ang_misset import get_item
    from LS49.ML_push.exploratory_missetting import metric
    R = get_item(idx)
    print ("coarse ground truth with index",idx)
    C = crystal_orientation.crystal_orientation(
      Amat,crystal_orientation.basis_type.direct)
    C.show(legend="ground truth, P1")
    C2 = C.change_basis(self.CB_OP_C_P.inverse())
    C2.show(legend="ground truth, C2")
    direct_A = R["integrated_crystal_model"].get_A_inverse_as_sqr() # from dials model, integration step
    permute = sqr((0,0,1,0,1,0,-1,0,0))
    sim_compatible = direct_A*permute # permute columns when post multiplying
    P = crystal_orientation.crystal_orientation(
      sim_compatible,crystal_orientation.basis_type.direct)
    P.show(legend="dials_integrated, C2")
    PR = P.change_basis(self.CB_OP_C_P)
    PR.show(legend="dials_integrated, primitive setting")
    PRC2 = PR.change_basis(self.CB_OP_C_P.inverse()) # dials integrated, C2
    cb_op_align = PR.best_similarity_transformation(C,200,1)
    align_PR = PR.change_basis(sqr(cb_op_align))
    align_PR.show(legend="dials_integrated, P1, aligned")
    print("alignment matrix", cb_op_align)
    metric_val = metric(align_PR,C)
    print("Key %d aligned angular offset is %12.9f deg."%(idx, metric_val))
    print("Key %d alignC2 angular offset is %12.9f deg."%(idx, metric(align_PR.change_basis(self.CB_OP_C_P.inverse()),C2)))
    from IPython import embed; embed()
    # coarse, dials crystal orientation models = C, align_PR
    # apply Rotx:
    import math
    align_PR_dx = align_PR.rotate_thru((1.0,0.0,0.0), math.pi* 0.01/180.)
    align_PR_dy = align_PR.rotate_thru((0.0,1.0,0.0), math.pi* 0.01/180.)
    align_PR_dz = align_PR.rotate_thru((0.0,0.0,1.0), math.pi* 0.01/180.)
    return (align_PR.direct_matrix(),align_PR_dx.direct_matrix(),
            align_PR_dy.direct_matrix(), align_PR_dz.direct_matrix()
"""
work out the basic code for
OK 1) getting CB_OP
OK 2) Getting sfall_main
OK 3) getting the coarse ground truth
OK 4) getting the dials refine
OK 5) applying rotational perturbations to dials refine
6) performing 7150 eV simulations for all three orientations (GT, RXGT, RYGT)
"""

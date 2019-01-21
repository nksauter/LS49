from __future__ import division, absolute_import, print_function
from scitbx.array_family import flex
#from dxtbx.format.Registry import Registry
from cctbx import crystal_orientation
import sys
import glob,math
import dials

# some parameters
json_glob = "/net/dials/raid1/sauter/LS49_integ_allrestr/idx*.img_integrated_experiments.json"
image_glob = "/net/dials/raid1/sauter/LS49/step5_MPIbatch_0%05d.img.gz"
global format_class

#integration_run = sys.argv[1] # LS49_integ_allrestr, LS49_integ_betarestr, LS49_integ_step5cori

#json_glob = json_glob.replace("LS49_integ_allrestr",integration_run)

# PART I. Generate read-header ABC matrices from img.gz files each run through
from LS49.work_pre_experiment.fine_detail_ground_truth import nanoBragg_mock
Mock_nano = nanoBragg_mock()
def get_items():
  file_list = glob.glob(json_glob)
  print ("There are %d items in the file list"%len(file_list))
  format_class = None
  for item in file_list:
    serial_no = int(item[-37:-32])
    image_file = image_glob%serial_no
    #print (image_file)
    if format_class is None:
      #format_class = Registry.find(image_file)
      from dxtbx.format.FormatSMVJHSim import FormatSMVJHSim
      format_class = FormatSMVJHSim
    i = format_class(image_file)
    Z = i.get_smv_header(image_file)
    ABC = Z[1]["DIRECT_SPACE_ABC"]
    abc = tuple([float(a) for a in ABC.split(",")])

    #from dxtbx.model.experiment_list import ExperimentListFactory
    #EC = ExperimentListFactory.from_json_file(item,check_format=False)[0].crystal

    print(abc[0],abc[1],abc[2])
    print(abc[3],abc[4],abc[5])
    print(abc[6],abc[7],abc[8])
    from LS49.work_pre_experiment.fine_detail_ground_truth import tst_all
    rotation = tst_all(serial_no)

    Mock_nano.set_rotation(rotation)
    FGabc = Mock_nano.get_average_abc()
    print(FGabc[0],FGabc[1],FGabc[2])
    print(FGabc[3],FGabc[4],FGabc[5])
    print(FGabc[6],FGabc[7],FGabc[8])
    yield(dict(serial_no=serial_no,fine=FGabc,coarse=abc))

#Part 2. Generate a look-up table of the ABC matrices
# dump these to pickle file:
#from six.moves import cPickle as pickle
#with open("dump_coarse_file.pickle","w") as out:
#  for stuff in get_items():
#    pickle.dump(stuff, out)
#exit("done")

#Part 3. Use a look-up table of the ABC matrices
def get_items():
  from six.moves import cPickle as pickle
  with open("dump_coarse_file.pickle","r") as inp:
    while 1:
      yield pickle.load(inp)

if __name__=="__main__":
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

    from cctbx import crystal_orientation

    # work up the ground truth from header
    header_Ori = crystal_orientation.crystal_orientation(stuff["coarse"], crystal_orientation.basis_type.direct)
    header_Ori.show(legend="header_Ori")

    C2_ground_truth = header_Ori.change_basis(CB_OP_C_P.inverse())
    C2_ground_truth.show(legend="coarse_ground_truth")

    # work up the ground truth from mosaic ensemble
    mosaic_Ori = crystal_orientation.crystal_orientation(stuff["fine"], crystal_orientation.basis_type.direct)
    mosaic_Ori.show(legend="mosaic_Ori")

    fine_ground_truth = mosaic_Ori.change_basis(CB_OP_C_P.inverse())
    fine_ground_truth.show(legend="fine_ground_truth")

    # now calculate the angle as mean a_to_a,b_to_b,c_to_c
    aoff = C2_ground_truth.a.angle(fine_ground_truth.a,deg=True)
    boff = C2_ground_truth.b.angle(fine_ground_truth.b,deg=True)
    coff = C2_ground_truth.c.angle(fine_ground_truth.c,deg=True)
    # solved:  the reason missetting_rot doesn't exactly align postref and ground_truth is
    # that it's a monoclinic lattice, not orthorhombic.  Difference in the beta angle prevents exact alignment

    hyp = flex.mean(flex.double((aoff,boff,coff)))

    angles.append(hyp)
    print("Item %5d serial_no %5d angular offset is %12.9f deg."%(icount,stuff["serial_no"],hyp))

  print("RMSD", math.sqrt(flex.mean(angles*angles)))

from __future__ import print_function, division
from six.moves import range
from six.moves import cPickle as pickle
from scitbx.array_family import flex
import scitbx
import math
from LS49.work2_for_aca_lsq.abc_background import fit_roi_multichannel # implicit import
# multichannel needed for unpickling

# %%% boilerplate specialize to packaged big data %%%
import os
from LS49.sim import step5_pad
from LS49.sim import step4_pad
from LS49.spectra import generate_spectra
ls49_big_data = os.environ["LS49_BIG_DATA"] # get absolute path from environment
step5_pad.big_data = ls49_big_data
step4_pad.big_data = ls49_big_data
generate_spectra.big_data = ls49_big_data
# %%%%%%
from LS49.sim.step5_pad import data
local_data = data()
Fe_oxidized_model = local_data.get("Fe_oxidized_model")
Fe_reduced_model = local_data.get("Fe_reduced_model")
Fe_metallic_model = local_data.get("Fe_metallic_model")
model_mode = os.environ["MODEL_MODE"]

def get_items(myrank,N_total,N_stride,cohort=0):
  #selected_idx = dict()
  #with open("trimmed_images") as F:
  #  for line in F:
  #    tokens = line.strip().split()
  #    angle = float(tokens[7]); idx = int(tokens[3])
  #    if angle < 0.028 : # just the slim angles
  #      selected_idx[idx]=angle
  for key in range(cohort*N_total, (cohort+1)*N_total):
    #each rank should only allow keys in the range from
    # cohort*N_total + myrank*N_stride to cohort*N_total + (myrank+1)*N_stride
    if key < cohort*N_total + myrank*N_stride: continue
    if key >= cohort*N_total + (myrank+1) * N_stride: continue
    if key >= (cohort+1) * N_total : continue
    #if key not in selected_idx:
    #  print("skipped image",key)
    #  continue
    try:
      with open("abc_coverage_%s/abcX%06d.pickle"%(model_mode,key),"rb") as F:
        T = pickle.load(F)
    except IOError:
      #print("No file abc_coverage/abcX%06d.pickle"%key)
      continue
    yield T,key


def pprint(M):
  islow,ifast=M.focus()
  for x in range(islow):
    print (" ".join([("%4.2f"%(10*M[(x,y)])) for y in range(ifast)]))



class MPI_Run(object):
  def __init__(self):
    from xfel.merging.application.mpi_helper import mpi_helper
    self.mpi_helper = mpi_helper()

  def __del__(self):
    self.mpi_helper.finalize()

  def parse_input(self):
    '''Parse input at rank 0 and broadcast the input parameters and options to all ranks'''

    if self.mpi_helper.rank == 0:
      from LS49.ML_push.phil import phil_scope
      help_message = '''Refine fp fdp parameters.'''

      # The script usage
      import libtbx.load_env
      self.usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name
      self.parser = None

      '''Initialize the script.'''
      from dials.util.options import OptionParser
      # Create the parser
      self.parser = OptionParser(
        usage=self.usage,
        phil=phil_scope,
        epilog=help_message)

      # Parse the command line. quick_parse is required for MPI compatibility
      params, options = self.parser.parse_args(show_diff_phil=True,quick_parse=True)

      # prepare for transmitting input parameters to all ranks
      transmitted = dict(params = params, options = options)
    else:
      transmitted = None

    # broadcast parameters and options to all ranks

    transmitted = self.mpi_helper.comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']

  def run(self):

    self.parse_input()

    N_total = self.params.N_total # number of items to simulate, nominally 100000
    logical_rank = self.mpi_helper.rank
    logical_size = self.mpi_helper.size
    if self.mpi_helper.rank==0 and self.mpi_helper.size==1: # special case of testing it
      try:
        logical_rank = self.params.tester.rank
        logical_size = self.params.tester.size
      except Exception: pass
    N_stride = int(math.ceil(N_total/logical_size)) # total number of tasks per rank
    print ("hello from rank %d of %d with stride %d"%(logical_rank,logical_size,N_stride))

    #from scitbx.lbfgs.tst_mpi_split_evaluator import mpi_split_evaluator_run
    #from scitbx.lbfgs.tst_mpi_split_evaluator import run_mpi as simple_tester
    #simple_tester()

    self.mpi_helper.comm.barrier()

    if self.mpi_helper.rank==0:
      print ("Finding initial G and abc factors")
    per_rank_items = []
    per_rank_keys = []
    per_rank_spots = 0
    per_rank_pixels = 0
    per_rank_millers = {}
    #per_rank_G = []
    min_spots = 3
    N_input=0
    from LS49.ML_push.new_global_fdp_refinery import fit_one_image_multispot
    for item,key in get_items(logical_rank,N_total,N_stride,self.params.cohort):
      N_input+=1
      if len(item) >= min_spots:
        #FOI = fit_one_image_multispot(list_of_images=item,
        #    HKL_lookup = transmitted_info["HKL_lookup"],
        #    model_intensities = transmitted_info["model_intensities"])

        #print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f"""%(
        #key, len(item), FOI.compute_functional_and_gradients()[0]))
        print ("""LLG Image %06d on %d Bragg spots NLL"""%(key, len(item)))
        # put the newly refined background model back into the item
        per_rank_items.append(item)
        per_rank_keys.append(key)

        #from IPython import embed; embed()
        nspots_this_image = len(per_rank_items[-1])
        per_rank_spots += nspots_this_image
        for ihkl in range(nspots_this_image):
          this_spot = per_rank_items[-1][ihkl]
          this_miller = this_spot.asu_idx_C2_setting
          #print(this_miller)
          per_rank_millers[this_miller]=1
          this_shoebox_size = len(this_spot.roi)
          per_rank_pixels += this_shoebox_size
        #  per_rank_items[-1][ihkl].bkgrd_a = flex.double(
        #            [FOI.a[3*ihkl+0],FOI.a[3*ihkl+1],FOI.a[3*ihkl+2]])
        #per_rank_G.append( FOI.a[-1] )

    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_items)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_items), self.mpi_helper.MPI.SUM, 0)
    N_input_images = self.mpi_helper.comm.reduce(N_input, self.mpi_helper.MPI.SUM, 0)
    N_shoeboxes = self.mpi_helper.comm.reduce(per_rank_spots, self.mpi_helper.MPI.SUM, 0)
    N_pixels = self.mpi_helper.comm.reduce(per_rank_pixels, self.mpi_helper.MPI.SUM, 0)
    all_miller_lists = self.mpi_helper.comm.gather(per_rank_millers, root=0)
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      print ("final report %d ranks, %d input images, %d refined images with >= 3 spots"%(
      N_ranks, N_input_images, N_refined_images))
      all_millers={}
      for item in all_miller_lists:
        all_millers.update(item)
      print ("final report %d shoeboxes, %d pixels, %d unique C2 asu Millers"%(
      N_shoeboxes, N_pixels, len(all_millers)))
      print ("Finished finding initial G and abc factors")
      print ("Initiating the full minimization")
    # -----------------------------------------------------------------------

    self.mpi_helper.comm.barrier()


if __name__=="__main__":
  Usage = """
mpirun -c 64 libtbx.python ../modules/LS49/ML_push/shoebox_stats.py N_total=8000
mpirun -c 64 libtbx.python ../modules/LS49/ML_push/shoebox_stats.py N_total=100000
  """

  script = MPI_Run()
  result = script.run()

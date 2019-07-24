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

from LS49.ML_push.new_global_fdp_refinery import get_items
from LS49.ML_push.differential_roi_manager import differential_roi_manager
from LS49.ML_push.shoebox_troubleshoot import pprint

#ersatz get_items for cori second run through to fill in gaps
def get_items(myrank,N_total,N_stride,shuffA,cohort=0):
  abc_glob = os.environ["ABC_GLOB"]
  abc_glob_pixel_ref = os.environ["ABC_GLOB_PIXEL_REF"]
  #shuffA = list(range(N_total))
  #import random
  #random.shuffle(shuffA)

  for key in range(cohort*N_total, (cohort+1)*N_total):
    #each rank should only allow keys in the range from
    # cohort*N_total + myrank*N_stride to cohort*N_total + (myrank+1)*N_stride
    if key < cohort*N_total + myrank*N_stride: continue
    if key >= cohort*N_total + (myrank+1) * N_stride: continue
    if key >= (cohort+1) * N_total : continue
    derkey = shuffA[key]
    abc_file_pixel_ref = abc_glob_pixel_ref%derkey
    if os.path.isfile(abc_file_pixel_ref): continue
    try:
      abc_file = abc_glob%derkey
      with open(abc_file,"rb") as F:
        T = pickle.load(F)
    except IOError:
      print("No file %s"%abc_file)
      continue
    yield T,derkey


abc_glob_pixel_ref = os.environ["ABC_GLOB_PIXEL_REF"] # output directory

class fit_one_image_multispot:
  def __init__(self,key,list_of_images,HKL_lookup,model_intensities,spectra,crystal):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    import scitbx
    #lay out the parameters.
    self.n_spots = len(list_of_images)
    self.n = 3*self.n_spots + 1 + 3
    self.x = flex.double()
    for ispot in range(self.n_spots):
      self.x.append(list_of_images[ispot].bkgrd_a[0])
      self.x.append(list_of_images[ispot].bkgrd_a[1])
      self.x.append(list_of_images[ispot].bkgrd_a[2])
    self.x.append(0.); self.x.append(0.); self.x.append(0.); # rotx, roty, rotz
    self.x.append(1.)
    # insert an ROI simulation here
    spectrum = spectra.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)
    self.DRM = differential_roi_manager(key,spotlist=list_of_images,spectrum=spectrum,crystal=crystal,allspectra = spectra)

    self.sb_data = []
    for ispot in range(self.n_spots):
      self.sb_data.append(list_of_images[ispot].sb_data)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=False,
        traditional_convergence_test_eps=0.05, #1.e-3 significantly (4x) quicker than 1.e-4
        drop_convergence_test_max_drop_eps=1.e-5, # using drop convergence test (since traditional = False), drastic effect
        #min_iterations=min_iterations,
        #max_iterations = None,
        max_calls=100)
    )
    self.a = self.x

  def update_roi_model_pixels_with_current_rotation(self,rotxyz):
    assert len(rotxyz) == 3
    updated_models4 = self.DRM.get_incremented_rotation_models(rotxyz)
    self.roi_model_pixels = []
    self.new_calc2_dict_last_round = []
    for ispot,spot in enumerate(self.list_of_images):

      M = self.DRM.data["miller_index"] # from dials integration pickle
      S = (spot.orig_idx_C2_setting) # happens to be the coarse ground truth spot
      idx = M.first_index(S)
      shoe = self.DRM.data["shoebox"][idx]
      B = shoe.bbox
      ROI = ((B[0],B[1]),(B[2],B[3]))
      new_calc2_dict = self.DRM.perform_one_simulation_optimized(model="Amat",ROI=ROI,models4 = updated_models4)

      intensity = new_calc2_dict["intensity"]
      this_P1_Miller_index = new_calc2_dict["miller"]
      #print("NEW",this_P1_Miller_index);pprint (new_calc2_dict["roi_pixels"])
      lookup_idx = self.HKL_lookup[this_P1_Miller_index]
      energy_dependent_intensity = self.model_intensities.matrix_copy_block(
        i_row=lookup_idx,i_column=0,n_rows=1,n_columns=100)
      rescale_factor = energy_dependent_intensity.as_1d() / intensity
      channels = new_calc2_dict["channels"]
      self.roi_model_pixels.append(rescale_factor[0] * channels[0])
      self.new_calc2_dict_last_round.append(new_calc2_dict) # need this after lbfgs when results are written out

      for ichannel in range(1,len(channels)):
        # rescale_factor[51] always == 1, equivalent to simtbx_intensity_7122
        self.roi_model_pixels[ispot] += rescale_factor[2*ichannel] * channels[2*ichannel]

      """
      intensity = list_of_images[ispot].simtbx_intensity_7122
      this_P1_Miller_index = list_of_images[ispot].simtbx_P1_miller
      print("2",this_P1_Miller_index)
      lookup_idx = HKL_lookup[this_P1_Miller_index]
      energy_dependent_intensity = model_intensities.matrix_copy_block(
        i_row=lookup_idx,i_column=0,n_rows=1,n_columns=100)
      rescale_factor = energy_dependent_intensity.as_1d() / intensity
      channels = list_of_images[ispot].channels
      self.roi_model_pixels[ispot]=rescale_factor[0] * channels[0]

      for ichannel in range(1,len(channels)):
        # rescale_factor[51] always == 1, equivalent to simtbx_intensity_7122
        self.roi_model_pixels[ispot] += rescale_factor[2*ichannel] * channels[2*ichannel]
      """
    # derivatives:
    self.roi_dxyz = dict(Amat_dx = [],Amat_dy = [], Amat_dz = [])
    for ispot,spot in enumerate(self.list_of_images):
      M = self.DRM.data["miller_index"] # from dials integration pickle
      S = (spot.orig_idx_C2_setting) # happens to be the coarse ground truth spot
      idx = M.first_index(S)
      shoe = self.DRM.data["shoebox"][idx]
      B = shoe.bbox
      ROI = ((B[0],B[1]),(B[2],B[3]))

      for modelkey in self.roi_dxyz:
        if modelkey == "Amat_dx": continue # don't consider rotX as it is parallel to the beam and well determined
        self.roi_dxyz[modelkey].append( -1. * self.roi_model_pixels[ispot] ) # finite diff derivative, subtract off base value

        new_dxyz_dict = self.DRM.perform_one_simulation_optimized(model=modelkey,ROI=ROI,models4 = updated_models4)

        intensity = new_dxyz_dict["intensity"]
        this_P1_Miller_index = new_dxyz_dict["miller"]
        lookup_idx = self.HKL_lookup[this_P1_Miller_index]
        energy_dependent_intensity = self.model_intensities.matrix_copy_block(
        i_row=lookup_idx,i_column=0,n_rows=1,n_columns=100)
        rescale_factor = energy_dependent_intensity.as_1d() / intensity
        channels = new_dxyz_dict["channels"]
        for ichannel in range(len(channels)):
          self.roi_dxyz[modelkey][ispot] += rescale_factor[2*ichannel] * channels[2*ichannel]

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%9.3f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
    self.update_roi_model_pixels_with_current_rotation(self.x[-4:-1]) # pass in rotational increments
    self.a = self.x
    f = 0.;
    g = flex.double(self.n)
    for ispot in range(self.n_spots):
      F = self.sb_data[ispot].focus()
      for x in range(F[1]):
        for y in range(F[2]):
          model_lambda = self.a[3*ispot+0]*x+self.a[3*ispot+1]*y+self.a[3*ispot+2]+ \
                         self.a[-1]*self.roi_model_pixels[ispot][x,y]
          datapt = self.sb_data[ispot][0,x,y] # not sure the right datapt when model_lambda<0
          if model_lambda<=0:
            f+= model_lambda # complete kludge, guard against math domain error
          else:
            f += model_lambda - datapt * math.log(model_lambda)
          one_minus_k_over_lambda = (1. - datapt/model_lambda)
          g[3*ispot+0] += x * one_minus_k_over_lambda # from handwritten notes
          g[3*ispot+1] += y * one_minus_k_over_lambda
          g[3*ispot+2] += one_minus_k_over_lambda
          # for this paper, never refine rotX as it is parallel to beam and well-determined
          g[-4] += one_minus_k_over_lambda * self.a[-1] * 0. # always fix rotx.  self.roi_dxyz["Amat_dx"][ispot][x,y]
          g[-3] += one_minus_k_over_lambda * self.a[-1] * self.roi_dxyz["Amat_dy"][ispot][x,y]
          g[-2] += one_minus_k_over_lambda * self.a[-1] * self.roi_dxyz["Amat_dz"][ispot][x,y]
          g[-1] += self.roi_model_pixels[ispot][x,y] * one_minus_k_over_lambda
    self.print_step("LBFGS stp",f)
    return f, g

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

    N_total = 100000 # self.params.N_total # number of items to simulate, nominally 100000
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

    if self.params.starting_model.algorithm=="to_file":
      if self.mpi_helper.rank == 0:
        HKL_lookup,static_fcalcs = get_static_fcalcs_with_HKL_lookup()
        from LS49.work2_for_aca_lsq.remake_range_intensities_with_complex \
           import get_intensity_structure
        model_intensities = get_intensity_structure(
           static_fcalcs,FE1_model=Fe_oxidized_model,FE2_model=Fe_reduced_model)
        with (open(self.params.starting_model.filename,"wb")) as out:
          pickle.dump(HKL_lookup,out,pickle.HIGHEST_PROTOCOL)
          pickle.dump(static_fcalcs,out,pickle.HIGHEST_PROTOCOL)
          pickle.dump(model_intensities,out,pickle.HIGHEST_PROTOCOL)
      return
    else:
      if self.mpi_helper.rank == 0:
        with (open(self.params.starting_model.filename,"rb")) as inp:
          print("the starting model (used for channel weighting) is",self.params.starting_model.filename)
          HKL_lookup = pickle.load(inp)
          static_fcalcs = pickle.load(inp)
          model_intensities = pickle.load(inp)
          from LS49.spectra.generate_spectra import spectra_simulation
          from LS49.sim.step5_pad import microcrystal

        shuffA = list(range(N_total))
        import random
        random.shuffle(shuffA)

        transmitted_info = dict(HKL_lookup = HKL_lookup,
                                static_fcalcs = static_fcalcs,
                                model_intensities = model_intensities,
                                spectra_simulation = spectra_simulation(),
                                crystal = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0),
                                shuffA = shuffA
                                )
      else:
        transmitted_info = None
    print ("before braodcast with ",self.mpi_helper.rank,self.mpi_helper.size)
    transmitted_info = self.mpi_helper.comm.bcast(transmitted_info, root = 0)
    self.mpi_helper.comm.barrier()
    print ("after barrier")

    # -----------------------------------------------------------------------
    if self.mpi_helper.rank==0:
      print ("Finding initial G and abc factors")
    per_rank_items = []
    per_rank_keys = []
    per_rank_G = []
    min_spots = 3
    N_input=0
    import os,omptbx # cori workaround, which does not get OMP_NUM_THREADS from environment
    workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
    omptbx.omp_set_num_threads(workaround_nt)
    for item,key in get_items(logical_rank,N_total,N_stride,transmitted_info["shuffA"],self.params.cohort):
      N_input+=1
      if len(item) >= min_spots:
        try:
          FOI = fit_one_image_multispot(key=key,
                                      list_of_images=item,
                                      HKL_lookup = transmitted_info["HKL_lookup"],
                                      model_intensities = transmitted_info["model_intensities"],
                                      spectra = transmitted_info["spectra_simulation"],
                                      crystal = transmitted_info["crystal"]
                                      )
        except RuntimeError as e:
          # no recovery from LBFGS error, skip event
          continue
        metric_P1, metric_C2 = FOI.DRM.get_current_angular_offsets(FOI.x[-4:-1])
        print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f angular offsets in P1 and C2 (degrees): %8.5f %8.5f"""%(
        key, len(item), FOI.compute_functional_and_gradients()[0],
        metric_P1, metric_C2),FOI.x[-4:-1]
        )

# reporting out results to new abc_coverage pickles.  Use the old one "item" as a template:
#    item := [<LS49.work2_for_aca_lsq.abc_background.fit_roi_multichannel>,... one for each spot]
# each fit_roi has the following attributes and we modify them as follows:
#        'a', the abcG parameters of the original one-spot fit, unused
#        'asu_idx_C2_setting', unmodified
#        'image_no', same as key, unmodified
#        'n', number of parameters for one-spot fit, unused
#        'orig_idx_C2_setting', unmodified
#        'sb_data', original shoebox data [integers as floats], passed along unmodified
#        'simtbx_P1_miller', unmodified
#        'simtbx_intensity_7122', unmodified
#        'x', the abcG parameters of the original one-spot fit, unused
# modify these:
#        'bkgrd_a', the spot abc parameters output here, to be passed on to global data fit
        for ispot in range(len(item)):
          item[ispot].bkgrd_a = FOI.x[3*ispot:3*(ispot+1)]
#        'channels', need to save the new data that was set during update_roi_model_pixels_with_current_rotation()
#                    but only for the Amat, not the derivatives.
          item[ispot].channels = FOI.new_calc2_dict_last_round[ispot]["channels"]
#        'roi', get new region of interest summation that was set during update_roi_model_pixels_with_current_rotation()
          item[ispot].roi = FOI.roi_model_pixels[ispot]

        # put the newly refined background model back into the item
        per_rank_keys.append(key)
        per_rank_G.append( FOI.a[-1] )

        print ("pickling modified abc_coverage file for key %d in rank %d"%(key,logical_rank),)
        with open(abc_glob_pixel_ref%(key),"wb") as F:
          pickle.dump(item,F, pickle.HIGHEST_PROTOCOL)

    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_keys)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_keys), self.mpi_helper.MPI.SUM, 0)
    N_input_images = self.mpi_helper.comm.reduce(N_input, self.mpi_helper.MPI.SUM, 0)
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      print ("final report %d ranks, %d input images, %d refined models"%(
      N_ranks, N_input_images, N_refined_images))
      print ("Finished finding initial G and abc factors")

    self.mpi_helper.comm.barrier()

if __name__=="__main__":
  Usage = """srun -n 32 -c 2 libtbx.python new_global_fdp_refinery.py #smallish test case, 1 node
mpirun -c 64 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py LLG_evaluator.enable_plot=True starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=8000

mpirun -c 56 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py LLG_evaluator.enable_plot=True starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=8000 LLG_evaluator.max_calls=19 LLG_evaluator.plot_interpolation=False

libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py LLG_evaluator.enable_plot=True # plot test
             ...either works only under: salloc -C haswell -N1 -q interactive -t 04:00:00
convert -delay 12 -loop 0 *.png Fe_metal_iteration.gif
  """

  script = MPI_Run()
  result = script.run()

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
from LS49.ML_push.shoebox_troubleshoot import pprint3,pprint

class fit_one_image_multispot:
  def __init__(self,key,list_of_images,HKL_lookup,model_intensities,spectra,crystal):
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
    self.roi_model_pixels = []
    # insert an ROI simulation here
    spectrum = spectra.generate_recast_renormalized_image(image=key,energy=7120.,total_flux=1e12)
    self.DRM = differential_roi_manager(key,spotlist=list_of_images,spectrum=spectrum,crystal=crystal,allspectra = spectra)
    #exit("for now, until I get ROI simulation going")

    for ispot in range(self.n_spots):
      intensity = list_of_images[ispot].simtbx_intensity_7122
      this_P1_Miller_index = list_of_images[ispot].simtbx_P1_miller
      lookup_idx = HKL_lookup[this_P1_Miller_index]
      energy_dependent_intensity = model_intensities.matrix_copy_block(
        i_row=lookup_idx,i_column=0,n_rows=1,n_columns=100)
      rescale_factor = energy_dependent_intensity.as_1d() / intensity
      channels = list_of_images[ispot].channels
      self.roi_model_pixels.append(rescale_factor[0] * channels[0])

      for ichannel in range(1,len(channels)):
        # rescale_factor[51] always == 1, equivalent to simtbx_intensity_7122
        self.roi_model_pixels[ispot] += rescale_factor[2*ichannel] * channels[2*ichannel]
      print ("in fit one image with",this_P1_Miller_index)
      pprint(self.roi_model_pixels[ispot])


here:
1) check that the printout here is the same as the differential_Roi_manager printout of dials_refine
2) check that the rescale_factor applied to live nanoBragg (differential_roi_manager) gives same output as printout here

    self.sb_data = []
    for ispot in range(self.n_spots):
      self.sb_data.append(list_of_images[ispot].sb_data)
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=0.05, #1.e-3 significantly (4x) quicker than 1.e-4
        #drop_convergence_test_max_drop_eps=max_drop_eps,
        #min_iterations=min_iterations,
        #max_iterations = None,
        max_calls=1000)
    )
    self.a = self.x

  def print_step(self,message,target):
    print ("%s %10.4f"%(message,target),
           "["," ".join(["%9.3f"%a for a in self.x]),"]")

  def compute_functional_and_gradients(self):
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
          g[3*ispot+0] += x * (1. - datapt/model_lambda) # from handwritten notes
          g[3*ispot+1] += y * (1. - datapt/model_lambda)
          g[3*ispot+2] += (1. - datapt/model_lambda)
          g[-1] += self.roi_model_pixels[ispot][x,y] * (1. - datapt/model_lambda)
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
          HKL_lookup = pickle.load(inp)
          static_fcalcs = pickle.load(inp)
          model_intensities = pickle.load(inp)
          from LS49.spectra.generate_spectra import spectra_simulation
          from LS49.sim.step5_pad import microcrystal

        transmitted_info = dict(HKL_lookup = HKL_lookup,
                                static_fcalcs = static_fcalcs,
                                model_intensities = model_intensities,
                                spectra_simulation = spectra_simulation(),
                                crystal = microcrystal(Deff_A = 4000, length_um = 4., beam_diameter_um = 1.0)
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
    for item,key in get_items(logical_rank,N_total,N_stride,self.params.cohort):
      N_input+=1
      if len(item) >= min_spots:
        FOI = fit_one_image_multispot(key=key,
                                      list_of_images=item,
                                      HKL_lookup = transmitted_info["HKL_lookup"],
                                      model_intensities = transmitted_info["model_intensities"],
                                      spectra = transmitted_info["spectra_simulation"],
                                      crystal = transmitted_info["crystal"]
                                      )

        print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f"""%(
        key, len(item), FOI.compute_functional_and_gradients()[0]))
        # put the newly refined background model back into the item
        per_rank_items.append(item)
        per_rank_keys.append(key)
        for ihkl in range(FOI.n_spots):
          per_rank_items[-1][ihkl].bkgrd_a = flex.double(
                    [FOI.a[3*ihkl+0],FOI.a[3*ihkl+1],FOI.a[3*ihkl+2]])
        per_rank_G.append( FOI.a[-1] )
      print("break for unit testing, only one image")
      break # for unit testing, only one image
    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_items)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_items), self.mpi_helper.MPI.SUM, 0)
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

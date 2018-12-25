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

from LS49.ML_push.new_global_fdp_refinery import MPI_Run as upstream_base_script
from LS49.ML_push.new_global_fdp_refinery import get_items, fit_one_image_multispot
from LS49.ML_push.new_global_fdp_refinery import rank_0_fit_all_f, george_sherrell_star

class MPI_Run(upstream_base_script):
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

    from scitbx.lbfgs.tst_mpi_split_evaluator import mpi_split_evaluator_run

    assert self.params.starting_model.algorithm!="to_file", \
           "run new_global__fdp_refinery.py first, to generate starting model"
    if self.mpi_helper.rank == 0:
        with (open(self.params.starting_model.filename,"rb")) as inp:
          HKL_lookup = pickle.load(inp)
          static_fcalcs = pickle.load(inp)
          model_intensities = pickle.load(inp)

        transmitted_info = dict(HKL_lookup = HKL_lookup,
        static_fcalcs = static_fcalcs, model_intensities = model_intensities)
    else:
        transmitted_info = None
    transmitted_info = self.mpi_helper.comm.bcast(transmitted_info, root = 0)
    self.mpi_helper.comm.barrier()
    # macrocycle 1 ---------------------------------------------------------
    # generate model_intensities table based on initial conditions
    if logical_rank == 0 or self.mpi_helper.size==1:
        FE1=local_data.get(self.params.starting_model.preset.FE1)
        FE2=local_data.get(self.params.starting_model.preset.FE2)

        from LS49.work2_for_aca_lsq.remake_range_intensities_with_complex \
           import get_intensity_structure
        new_model_intensities = get_intensity_structure(
           transmitted_info["static_fcalcs"],FE1_model=FE1,FE2_model=FE2)
        broadcast_info = new_model_intensities
    else:
        broadcast_info = None
    current_model_intensities = self.mpi_helper.comm.bcast(broadcast_info, root = 0)
    self.mpi_helper.comm.barrier()


    # -----------------------------------------------------------------------
    if self.mpi_helper.rank==0:
      print ("Finding initial G and abc factors")
    per_rank_items = []
    per_rank_keys = []
    per_rank_G = []
    min_spots = 3
    N_input=0
    for item,key in get_items(logical_rank,N_total,N_stride):
      N_input+=1
      if len(item) >= min_spots:
        per_rank_items.append(item)
        per_rank_keys.append(key)
        FOI = fit_one_image_multispot(list_of_images=item,
            HKL_lookup = transmitted_info["HKL_lookup"],
            model_intensities = current_model_intensities)

        print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f"""%(
        key, len(item), FOI.compute_functional_and_gradients()[0]))
        # put the newly refined background model back into the item
        for ihkl in range(FOI.n_spots):
          per_rank_items[-1][ihkl].bkgrd_a = flex.double(
                    [FOI.a[3*ihkl+0],FOI.a[3*ihkl+1],FOI.a[3*ihkl+2]])
        per_rank_G.append( FOI.a[-1] )

    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_items)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_items), self.mpi_helper.MPI.SUM, 0)
    N_input_images = self.mpi_helper.comm.reduce(N_input, self.mpi_helper.MPI.SUM, 0)
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      print ("final report %d ranks, %d input images, %d refined models"%(
      N_ranks, N_input_images, N_refined_images))
      print ("Finished finding initial G and abc factors")
      print ("Initiating the full minimization")
    # -----------------------------------------------------------------------

    W = rank_0_fit_all_f( self.params,
                          FE1_model=local_data.get(self.params.starting_model.preset.FE1),
                          FE2_model=local_data.get(self.params.starting_model.preset.FE2))
    W.reinitialize(logical_rank, self.mpi_helper.size, per_rank_items, per_rank_keys, per_rank_G,
                 transmitted_info["HKL_lookup"],
                 transmitted_info["static_fcalcs"],current_model_intensities,
                 force_recompute=True)
    W.set_macrocycle(1)
    minimizer = mpi_split_evaluator_run(target_evaluator=W,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-2,
        max_calls=self.params.LLG_evaluator.max_calls)
      )
    if logical_rank==0:
      print("Minimizer ended at iteration",W.iteration)
    self.mpi_helper.comm.barrier()
    # 2nd macrocycle---------------------------------------------------------
    # generate model_intensities table based on initial conditions
    FE1 = george_sherrell_star(fp = W.x[0:100],fdp = W.x[100:200])
    FE2 = george_sherrell_star(fp = W.x[200:300],fdp = W.x[300:400])
    if logical_rank == 0 or self.mpi_helper.size==1:

        from LS49.work2_for_aca_lsq.remake_range_intensities_with_complex \
           import get_intensity_structure
        new_model_intensities = get_intensity_structure(
           transmitted_info["static_fcalcs"],FE1_model=FE1,FE2_model=FE2)
        broadcast_info = new_model_intensities
    else:
        broadcast_info = None
    current_model_intensities = self.mpi_helper.comm.bcast(broadcast_info, root = 0)
    self.mpi_helper.comm.barrier()


    # -----------------------------------------------------------------------
    if self.mpi_helper.rank==0:
      print ("Finding initial G and abc factors")
    per_rank_items = []
    per_rank_keys = []
    per_rank_G = []
    min_spots = 3
    N_input=0
    for item,key in get_items(logical_rank,N_total,N_stride):
      N_input+=1
      if len(item) >= min_spots:
        per_rank_items.append(item)
        per_rank_keys.append(key)
        FOI = fit_one_image_multispot(list_of_images=item,
            HKL_lookup = transmitted_info["HKL_lookup"],
            model_intensities = current_model_intensities)

        print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f"""%(
        key, len(item), FOI.compute_functional_and_gradients()[0]))
        # put the newly refined background model back into the item
        for ihkl in range(FOI.n_spots):
          per_rank_items[-1][ihkl].bkgrd_a = flex.double(
                    [FOI.a[3*ihkl+0],FOI.a[3*ihkl+1],FOI.a[3*ihkl+2]])
        per_rank_G.append( FOI.a[-1] )

    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_items)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_items), self.mpi_helper.MPI.SUM, 0)
    N_input_images = self.mpi_helper.comm.reduce(N_input, self.mpi_helper.MPI.SUM, 0)
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      print ("final report %d ranks, %d input images, %d refined models"%(
      N_ranks, N_input_images, N_refined_images))
      print ("Finished finding initial G and abc factors")
      print ("Initiating the full minimization")
    # -----------------------------------------------------------------------
    W_previous = W
    W = rank_0_fit_all_f( self.params,FE1_model=FE1,FE2_model=FE2)
    W.reinitialize(logical_rank, self.mpi_helper.size, per_rank_items, per_rank_keys, per_rank_G,
                 transmitted_info["HKL_lookup"],
                 transmitted_info["static_fcalcs"],current_model_intensities,
                 force_recompute=True)
    W.set_macrocycle(2, W_previous.starting_params_FE1, W_previous.starting_params_FE2)
    minimizer = mpi_split_evaluator_run(target_evaluator=W,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-2,
        max_calls=self.params.LLG_evaluator.max_calls)
      )
    if logical_rank==0:
      print("Minimizer ended at iteration",W.iteration)
    self.mpi_helper.comm.barrier()
    # 3rd macrocycle-----------------------------------------------------------
    # generate model_intensities table based on initial conditions
    FE1 = george_sherrell_star(fp = W.x[0:100],fdp = W.x[100:200])
    FE2 = george_sherrell_star(fp = W.x[200:300],fdp = W.x[300:400])
    if logical_rank == 0 or self.mpi_helper.size==1:

        from LS49.work2_for_aca_lsq.remake_range_intensities_with_complex \
           import get_intensity_structure
        new_model_intensities = get_intensity_structure(
           transmitted_info["static_fcalcs"],FE1_model=FE1,FE2_model=FE2)
        broadcast_info = new_model_intensities
    else:
        broadcast_info = None
    current_model_intensities = self.mpi_helper.comm.bcast(broadcast_info, root = 0)
    self.mpi_helper.comm.barrier()


    # -----------------------------------------------------------------------
    if self.mpi_helper.rank==0:
      print ("Finding initial G and abc factors")
    per_rank_items = []
    per_rank_keys = []
    per_rank_G = []
    min_spots = 3
    N_input=0
    for item,key in get_items(logical_rank,N_total,N_stride):
      N_input+=1
      if len(item) >= min_spots:
        per_rank_items.append(item)
        per_rank_keys.append(key)
        FOI = fit_one_image_multispot(list_of_images=item,
            HKL_lookup = transmitted_info["HKL_lookup"],
            model_intensities = current_model_intensities)

        print ("""LLG Image %06d on %d Bragg spots NLL    channels F = %9.1f"""%(
        key, len(item), FOI.compute_functional_and_gradients()[0]))
        # put the newly refined background model back into the item
        for ihkl in range(FOI.n_spots):
          per_rank_items[-1][ihkl].bkgrd_a = flex.double(
                    [FOI.a[3*ihkl+0],FOI.a[3*ihkl+1],FOI.a[3*ihkl+2]])
        per_rank_G.append( FOI.a[-1] )

    print ("rank %d has %d refined images"%(logical_rank,len(per_rank_items)))

    N_ranks = self.mpi_helper.comm.reduce(1, self.mpi_helper.MPI.SUM, 0)
    N_refined_images = self.mpi_helper.comm.reduce(len(per_rank_items), self.mpi_helper.MPI.SUM, 0)
    N_input_images = self.mpi_helper.comm.reduce(N_input, self.mpi_helper.MPI.SUM, 0)
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      print ("final report %d ranks, %d input images, %d refined models"%(
      N_ranks, N_input_images, N_refined_images))
      print ("Finished finding initial G and abc factors")
      print ("Initiating the full minimization")
    # -----------------------------------------------------------------------
    W_previous = W
    W = rank_0_fit_all_f( self.params,FE1_model=FE1,FE2_model=FE2)
    W.reinitialize(logical_rank, self.mpi_helper.size, per_rank_items, per_rank_keys, per_rank_G,
                 transmitted_info["HKL_lookup"],
                 transmitted_info["static_fcalcs"],current_model_intensities,
                 force_recompute=True)
    W.set_macrocycle(3, W_previous.starting_params_FE1, W_previous.starting_params_FE2)
    minimizer = mpi_split_evaluator_run(target_evaluator=W,
        termination_params=scitbx.lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-2,
        max_calls=self.params.LLG_evaluator.max_calls)
      )
    if logical_rank==0:
      print("Minimizer ended at iteration",W.iteration)
    self.mpi_helper.comm.barrier()

  """
1.5) When using wrong-redox initial state, not right to use pre-packaged model intensities, rather on-the-fly
3) use macrocycle over abcG / fdp refinement
Scale up to 12000 images, 64 cores, iteration to convergence
starting points with all four pre-set ox/red states
set regular list of trials needed for publication
Migrate to cori
  """
if __name__=="__main__":
  Usage = """srun -n 32 -c 2 libtbx.python macrocycle_refinery.py #smallish test case, 1 node
mpirun -c 64 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py LLG_evaluator.enable_plot=True starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=8000

mpirun -c 56 libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py LLG_evaluator.enable_plot=True starting_model.preset.FE1=Fe_metallic_model starting_model.preset.FE2=Fe_metallic_model N_total=2000 LLG_evaluator.max_calls=12 LLG_evaluator.plot_interpolation=False LLG_evaluator.title=macrocycle_2000

libtbx.python ../modules/LS49/ML_push/macrocycle_refinery.py LLG_evaluator.enable_plot=True # plot test
             ...either works only under: salloc -C haswell -N1 -q interactive -t 04:00:00
convert -delay 12 -loop 0 *.png Fe_metal_iteration.gif
  """

  script = MPI_Run()
  result = script.run()

Dec. 21, 2018

A calculation done to implement Max Likelihood refinement of the fp, fdp parameters.

Step 1: in working directory, remake the confirm_sfall_P1_7122_amplitudes.pickle
  libtbx.python ../modules/LS49/work2_for_aca_lsq/remake_7122_intensities.py

Step 2: tests
  libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py 45 # tests rank 45 
  mpirun -c 32 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py # actual MPI use

  # calculation of the first-time intensities, long way around
  mpirun -c 32 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py \ 
    starting_model.algorithm=to_file  --> output file new_global_fdp_big_data.pickle


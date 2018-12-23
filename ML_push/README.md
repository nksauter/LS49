Nov. 23

Looks like the NERSC scratch obliterated a lot of my July 4 result files.  

Pop the stack:
can't run refinery.py without re-generating abc_background.py
can't do that without re-running remake_7122_intensities.py

Refactor all of these scripts with the ls49_big_data model

Can't grok the plot.  Refactor sim/fdp_plot.py
libtbx.python ../modules/LS49/work_pre_experiment/interpolate_fdp_plot.py # generate 1-eV spacings for metal Fe

edit Fe-fake.dat

libtbx.python ../modules/LS49/work2_for_aca_lsq/remake_7122_intensities.py

libtbx.python ../modules/LS49/work2_for_aca_lsq/remake_range_intensities.py

libtbx.python ../modules/LS49/work2_for_aca_lsq/abc_background.py

Had to restore integration results and merging results from dials.lbl.gov:
scp -pr sauter@dials.lbl.gov:/net/dials/raid1/sauter/LS49_merge .
rsync -r sauter@dials.lbl.gov:/net/dials/raid1/sauter/LS49_integ_step5cori ./
(Then linked it to the LS49_integ_step5 directory)

finally refinery.py followed by plotLLG.py to reproduce July 4-5 results.
-----------------------------------------------------
Dec. 21, 2018.  Move over to dials.lbl.gov:/dev/shm
rsync -r nksauter@cori.nersc.gov:/global/cscratch1/sd/nksauter/proj-h0918/ls49nov/abc_coverage ./

A calculation done to implement Max Likelihood refinement of the fp, fdp parameters.

Step 1: in working directory, remake the confirm_sfall_P1_7122_amplitudes.pickle
  libtbx.python ../modules/LS49/work2_for_aca_lsq/remake_7122_intensities.py

Step 2: tests
  libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py 45 # tests rank 45 
  mpirun -c 32 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py # actual MPI use

  # calculation of the first-time intensities, long way around
  mpirun -c 32 libtbx.python ../modules/LS49/ML_push/new_global_fdp_refinery.py \ 
    starting_model.algorithm=to_file  --> output file new_global_fdp_big_data.pickle

bash
export WORK=/dev/shm/sauter
cd ${WORK}
source miniconda2/etc/profile.d/conda.sh
conda activate myEnv
export LS49_BIG_DATA=${WORK}/ls49_big_data
export OMP_NUM_THREADS=64
source build/setpaths.sh

Installation:
used junit-xml=1.7 instead of 1.8
conda install mpich2
conda install mpi openmpi mpi4py
conda install mrcfile
mpirun -c 32 libtbx.python modules/cctbx_project/scitbx/lbfgs/tst_mpi_split_eval
uator.py

cd ${WORK}/test
rm -rf ${WORK}/test/*
libtbx.python ${WORK}/modules/LS49/tests/public-test-all.py
rm -rf ${WORK}/test/*
libtbx.run_tests_parallel module=LS49 nproc=12
git branch ml_push
git push -u origin ml_push


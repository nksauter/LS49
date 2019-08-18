
```
export WORK=/global/cscratch1/sd/nksauter/proj-2col
cd $WORK
cp ../ls49.tgz.bolotovsky ./ls49.tgz
tar xzvf ls49.tgz
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate
chmod +x bootstrap.py
./bootstrap.py hot update --builder=dials
cd ${WORK}/modules
git clone git@github.com:nksauter/LS49.git
cd ${WORK}
```

```
module purge; module load python esslurm
salloc -C gpu -N 1 -A m1759 -t 04:00:00 --gres=gpu:1 -c 10
module load python/2.7-anaconda-2019.07
module swap PrgEnv-intel PrgEnv-gnu
conda create -n proj2col --clone base

module load cuda

source /usr/common/software/python/2.7-anaconda-2019.07/etc/profile.d/conda.sh
conda activate proj2col
pip install procrunner tqdm orderedset
conda install -c bkpoon wxpython

wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
pip uninstall mpi4py
tar zxf mpi4py-3.0.0.tar.gz
cd mpi4py-3.0.0/
module load gcc/7.3.0 cuda mvapich2 # or pgi/intel instead of gcc)
python setup.py build --mpicc=mpicc
python setup.py install
cd ${WORK}; srun -n 5 -c 2 python test.py # simple test of COMM_WORLD and get_rank()

mkdir ${WORK}/build
cd ${WORK}/build
python ../modules/cctbx_project/libtbx/configure.py LS49 prime iota --enable_cuda --enable_openmp_if_possible=True --use_conda

source setpaths.sh
make # or libtbx.scons -j 10 (repeated twice)

export OMP_NUM_THREADS=24
export LS49_BIG_DATA=${WORK}/ls49_big_data 
mkdir ${WORK}/xxx; cd ${WORK}/xxx
# comment out the tests on tst_polychromatic_image and double_precision_poly, as the OpenMP seems not to work on Cori
srun libtbx.run_tests_parallel module=LS49 nproc=13 # test to make sure regression test runs

# Adapt the step5_batch.py test
# set N_total to be the total number of images to simulate
# in step5_pad.py change the add_spots_algorithm to "cuda"
# other new edits to the code to specify device_id
# not completely checked in yet

### Each time login
module purge; module load python esslurm
salloc -C gpu -N 1 -A m1759 -t 04:00:00 --gres=gpu:1 -c 10
source /usr/common/software/python/2.7-anaconda-2019.07/etc/profile.d/conda.sh
conda activate proj2col
export WORK=/global/cscratch1/sd/nksauter/proj-2col
cd $WORK
export LS49_BIG_DATA=${WORK}/ls49_big_data
export OMP_NUM_THREADS=24
source build/setpaths.sh
module load gcc/7.3.0 cuda mvapich2
```


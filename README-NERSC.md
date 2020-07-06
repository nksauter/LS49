## Quick start instructions for building and executing, on NERSC platforms

### General setup (for all three platforms):  
Move ~/.condarc out of path  
create a working directory ${WORK}  
obtain the big data file ls49.tgz from nksauter
```
cd ${WORK}
tar zxf ls49.tgz
export LS49_BIG_DATA=${WORK}/ls49_big_data 
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate
python bootstrap.py hot update --builder=dials
cd ${WORK}/modules
git clone git@github.com:nksauter/LS49.git
cd ${WORK}
wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
```

### Software build on EDISON
```
bash
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
conda create -n EDISONapr24 --clone base
source /global/common/edison/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate EDISONapr24
pip install future mrcfile procrunner tqdm orderedset
conda install wxpython=3
pip uninstall mpi4py
tar zxf mpi4py-3.0.0.tar.gz
cd mpi4py-3.0.0/
LDFLAGS="-shared" python setup.py build --mpicc=$(which cc)
python setup.py build_exe --mpicc="$(which cc) -dynamic"
python setup.py install
python setup.py install_exe
cd ..
mkdir ${WORK}/build_edison_apr24
cd ${WORK}/build_edison_apr24
python ../modules/cctbx_project/libtbx/configure.py LS49 prime iota --enable_openmp_if_possible=True --use_conda
source setpaths.sh
make
export OMP_NUM_THREADS=24
mkdir ${WORK}/xxx; cd xxx
libtbx.run_tests_parallel module=LS49 nproc=20 # test to make sure regression test runs
```

### Subsequent login on EDISON:
```
salloc -A m2859 -C -N1 -q debug -t 00:30:00
bash
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/edison/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate EDISONapr24
source build_edison_apr24/setpaths.sh
export OMP_NUM_THREADS=24
export LS49_BIG_DATA=/global/cscratch1/sd/nksauter/proj-h0918/ls49_big_data 
```



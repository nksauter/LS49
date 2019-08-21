## Instructions for building and executing the benchmark:

### General setup (for all three platforms):  
Move ~/.condarc out of path  
create a working directory ${WORK}  
obtain the big data file ls49.tgz
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

### Exact software version:
Once source update was complete, the exact version for each package dependency was determined with "git log" and "git remote show origin":

| package        | remote origin           | guid                                          |
|----------------|----------------------------------|---------------------------------------------------|
| LS49           | git@github.com:nksauter/LS49.git      | 08f11d4ea399a882c7137f974aade58696e984e6 | 
| annlib_adaptbx | https://github.com/cctbx/annlib_adaptbx.git | fda57b44c32927313ba322de2ad3ff6d77a722fd |
| boost          | https://github.com/cctbx/boost.git | 72edfa007c73dd3dd01f5530205f7ca14cda7420 |
| cbflib         | https://github.com/yayahjb/cbflib.git | 611f21bbef81c5104954aaab9b72e4e55ec3e70d |
| cctbx_project  | https://github.com/cctbx/cctbx_project.git | 130a6cd07c60b563e656c339d52b3907c429ed81 |
| dials          | https://github.com/dials/dials.git| 1a873c105486a8e5923c82c540c6c9c848300bb0 |
| dxtbx          | https://github.com/cctbx/dxtbx.git | 8094491e1b239c9c16543ea7e0e1b47e40b2c8a5|
| xia2           | https://github.com/xia2/xia2.git | 082c6442f56cecc1ccd2933d169ee47a90aa4504|


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
mkdir ${WORK}/xxx; cd ${WORK}/xxx
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
export LS49_BIG_DATA=${WORK}/ls49_big_data 
```

### Run the benchmark on EDISON:
```
mkdir ${WORK}/050; cd ${WORK}/050
sbatch slurm_edison.sh
```

### Python environment build on Cori, good for both Haswell and KNL
```
salloc -A m2859 -C knl -N1 -q interactive -t 04:00:00
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/cori/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda create -n KNLapr22 --clone base
conda activate KNLapr22
pip install future mrcfile procrunner tqdm orderedset
conda install wxpython=3
pip uninstall mpi4py
tar zxf mpi4py-3.0.0.tar.gz
cd mpi4py-3.0.0/
python setup.py build --mpicc=$(which cc)
python setup.py build_exe --mpicc="$(which cc) -dynamic"
python setup.py install
python setup.py install_exe
cd ..
```

### Software build on Haswell
```
salloc -A m2859 -C knl -N1 -q interactive -t 04:00:00
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/cori/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate KNLapr22
mkdir build_haswell_may02
cd build_haswell_may02
python ../modules/cctbx_project/libtbx/configure.py LS49 prime iota --enable_openmp_if_possible=True --use_conda
source setpaths.sh
make
export OMP_NUM_THREADS=24
export LS49_BIG_DATA=${WORK}/ls49_big_data 
mkdir ${WORK}/xxx; cd ${WORK}/xxx
libtbx.run_tests_parallel module=LS49 nproc=20 # test to make sure regression test runs
```

### Run the benchmark on cori/Haswell:
```
salloc -A m2859 -C haswell -N1 -q interactive -t 04:00:00
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/cori/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate KNLapr22source build_edison_apr24/setpaths.sh
export OMP_NUM_THREADS=24
export LS49_BIG_DATA=${WORK}/ls49_big_data 
source build_haswell_may02/setpaths.sh
mkdir ${WORK}/051; cd ${WORK}/051
sbatch slurm_haswell.sh
```

### Software build on KNL
```
salloc -A m2859 -C knl -N1 -q interactive -t 04:00:00
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/cori/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate KNLapr22
mkdir build_knl_apr24
cd build_knl_apr24
python ../modules/cctbx_project/libtbx/configure.py LS49 prime iota --enable_openmp_if_possible=True --use_conda
source setpaths.sh
make
export OMP_NUM_THREADS=24
export LS49_BIG_DATA=${WORK}/ls49_big_data 
mkdir ${WORK}/xxx; cd ${WORK}/xxx
libtbx.run_tests_parallel module=LS49 nproc=20 # test to make sure regression test runs
```

### Run the benchmark on cori/KNL:
```
salloc -A m2859 -C knl -N1 -q interactive -t 04:00:00
cd ${WORK}
module load python/2.7-anaconda-5.2
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
source /global/common/cori/software/python/2.7-anaconda-5.2/etc/profile.d/conda.sh
conda activate KNLapr22
source build_knl_apr24/setpaths.sh
export OMP_NUM_THREADS=24
export LS49_BIG_DATA=${WORK}/ls49_big_data 
mkdir ${WORK}/052; cd ${WORK}/052
sbatch slurm_knl.sh
```

### Notes on connecting to Gitlab
```
eval `ssh-agent` # start agent, git port id, alternatively ssh-agent bash
ssh-add ~/.ssh/id_rsa # name the key, give passphrase
git clone git@gitlab.com:NESAP/exafel/cctbx/benchmark.git
```


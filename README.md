## Quick start instructions for building and executing, on linux command line

### General setup:  
Move ~/.condarc out of path  
create a working directory ${WORK}  
```
cd ${WORK}
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py --no-check-certificate
python bootstrap.py --builder=dials --mpi-build --use-conda --python=37 \
--config_flags="--enable_openmp_if_possible=True" --nproc=64 # conform to the number of cores available
cd ${WORK}/modules
git clone git@github.com:nksauter/LS49.git
git clone https://gitlab.com/cctbx/ls49_big_data.git

source ${WORK}/build/conda_setpaths.sh
conda install git git-lfs
cd ${WORK}/modules/ls49_big_data
git lfs install --local
git lfs pull  # gets the big data
source ${WORK}/build/conda_unsetpaths.sh
``` 

### Software build for linux
```
cd ${WORK}/build
source setpaths.sh
libtbx.configure LS49
make
export OMP_NUM_THREADS=24
mkdir ${WORK}/xxx; cd xxx
libtbx.run_tests_parallel module=LS49 nproc=20 # test to make sure regression test runs
```

### Subsequent login for linux
```
cd ${WORK}
source build/setpaths.sh
export OMP_NUM_THREADS=24
```

### Testing CUDA on Linux 
```
# at the bootstrap step (above) add cuda: --config_flags="--enable_cuda"
# On Cori-GPU, run the tests inside of an srun command (allocates GPU):
srun -n 1 -c 10 libtbx.run_tests_parallel nproc=Auto module=LS49 module=simtbx
```


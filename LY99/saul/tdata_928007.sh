#!/bin/bash -l
#SBATCH -N 4               # Number of nodes
#SBATCH -J tdata
#SBATCH -L SCRATCH          # job requires SCRATCH files
#SBATCH -A m3890_g          # allocation
#SBATCH -C gpu
#SBATCH -q early_science    # regular queue
#SBATCH -t 00:20:00         # wall clock time limit
#SBATCH -o job%j.out
#SBATCH -e job%j.err

export WORK=$SCRATCH/adse13_249/LY99
cd $WORK

mkdir $SLURM_JOB_ID; cd $SLURM_JOB_ID

export TRIAL=ly99sim
export OUT_DIR=${PWD}
# NO PSF:
export DIALS_OUTPUT=${WORK}/927185
# WITH PSF:
#export DIALS_OUTPUT=${WORK}/927187

echo "dispatch.step_list=input tdata
input.path=${DIALS_OUTPUT}
input.experiments_suffix=_integrated.expt
input.reflections_suffix=_integrated.refl
input.parallel_file_load.method=uniform
tdata.output_path=${TRIAL}_cells
output.prefix=${TRIAL}
output.output_dir=${OUT_DIR}/${TRIAL}/out
output.tmp_dir=${OUT_DIR}/${TRIAL}/tmp
output.do_timing=True
output.log_level=0
" > tdata.phil

mkdir -p ${OUT_DIR}/${TRIAL}/out
mkdir -p ${OUT_DIR}/${TRIAL}/tmp

echo "jobstart $(date)";pwd
srun -n 256 -c 2 cctbx.xfel.merge tdata.phil
echo "jobend $(date)";pwd

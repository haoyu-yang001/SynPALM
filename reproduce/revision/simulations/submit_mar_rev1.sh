#!/bin/bash
# S5 follow-up (rep_sims_mar_rev1.R): ./submit_mar_rev1.sh <S5W|S5V> <setting> [first] [last] [n_rep]
#SBATCH --job-name=rev1_mar
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/mar_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/mar_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=0-06:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --job-name=rev1_${1}_s$2 --array=${3:-1}-${4:-20} --export=ALL,STUDY=$1,SET=$2,NREP=${5:-25} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/rep_sims_mar_rev1.R $STUDY $SET $SLURM_ARRAY_TASK_ID $NREP

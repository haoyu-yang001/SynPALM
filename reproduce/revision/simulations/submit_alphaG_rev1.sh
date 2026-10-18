#!/bin/bash
# S2 alpha_G-only SNPs: ./submit_alphaG_rev1.sh <setting> [first] [last] [n_rep]
#SBATCH --job-name=rev1_s2
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/s2_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/s2_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=2-00:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --job-name=rev1_s2_s$1 --array=${2:-1}-${3:-100} --export=ALL,SET=$1,NREP=${4:-20} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/alphaG_rev1.R $SET $SLURM_ARRAY_TASK_ID $NREP

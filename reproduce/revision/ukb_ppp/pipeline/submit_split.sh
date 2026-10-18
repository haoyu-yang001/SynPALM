#!/bin/bash
# Split-sample replication (no PRS/PCs) for proteins on lines <start>-<end> of
# order_split_top100.txt (top 100 by out-of-fold rho).  ./submit_split.sh 1 100 [max_parallel]
#SBATCH --job-name=cf3_split
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_split/sp_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_split/sp_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=0-10:00
BASE=<PROJECT_DIR>
ORDER_FILE=${ORDER_FILE:-$BASE/code/sh_run/run_noPRS_noPC/order_split_top100.txt}
mkdir -p $BASE/logs/noPRS_noPC_split
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=${1:-1}-${2:-100}%${3:-100} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
p=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ORDER_FILE" | tr -d '[:space:]')
Rscript code/synpalm_main_code/split_validation_crossfit.R $p 10

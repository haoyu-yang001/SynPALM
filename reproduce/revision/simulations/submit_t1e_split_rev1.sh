#!/bin/bash
# Type I error with relatives split between labelled and unlabelled sets (POP-GWAS assumption violated):
#   ./submit_t1e_split_rev1.sh [first=1] [last=200] [n_chunks=250] [miss list, colon-separated, default 0.5:0.75]
#SBATCH --job-name=rev1_t1esplit
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/t1esplit_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/t1esplit_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2-00:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=${1:-1}-${2:-200} --export=ALL,NCH=${3:-250},MISS=${4:-0.5:0.75} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/typeIerror_split_rev1.R $SLURM_ARRAY_TASK_ID $NCH ${MISS//:/,}

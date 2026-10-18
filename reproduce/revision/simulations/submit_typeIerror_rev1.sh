#!/bin/bash
# Type I error simulations (revision): ./submit_typeIerror_rev1.sh <GRMr> [first] [last] [n_chunks] [miss list, colon-separated e.g. 0:0.5:0.75]
#   GRMr 0.5 (siblings, Table 1 incl. POP-GWAS, Fig 2B, Fig S4), 0.25 (half-siblings), 0.125 (cousins)
#SBATCH --job-name=rev1_t1e
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/t1e_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/t1e_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=3-00:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --job-name=rev1_t1e_$1_${5:-all} --array=${2:-1}-${3:-500} --export=ALL,GRMR=$1,NCH=${4:-1000},MISS=${5:-0:0.25:0.5:0.75:0.9} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/typeIerror_rev1.R $GRMR $SLURM_ARRAY_TASK_ID $NCH ${MISS//:/,}

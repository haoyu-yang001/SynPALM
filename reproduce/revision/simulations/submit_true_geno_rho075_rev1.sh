#!/bin/bash
# rho = 0.75 rerun of the true-genotype (UKB) simulations:
#   (t1e mode SUPERSEDED by submit_true_geno_t1e_calgeno_rev1.sh; not used)
#   ./submit_true_geno_rho075_rev1.sh t1e   [first] [last] [n_arrays=1000] [miss list, colon-separated]
#       true_geno_typeIerror_rho075_rev1.R: n_arrays tasks, each ~575 autosomal SNPs x 5 missing rates
#   ./submit_true_geno_rho075_rev1.sh power [first] [last] [miss list, colon-separated] [rho=0.75]
#       true_geno_power_rho075_rev1.R: 1000 tasks (20 sampled SNPs each) x 5 missing rates x 7 pve_g
# Tasks skip outputs that already exist, so requeued / resubmitted tasks resume.
#SBATCH --job-name=rev1_tg
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/tg_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/tg_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-12:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  MODE=$1
  if [ "$MODE" = t1e ]; then
    sbatch --job-name=rev1_tg_t1e_${5:-all} --mem=12G --time=0-04:00 --array=${2:-1}-${3:-1000} \
      --export=ALL,MODE=t1e,NARR=${4:-1000},MISS=${5:-0:0.25:0.5:0.75:0.9} "$0"
  else
    sbatch --job-name=rev1_tg_power_r${5:-0.75} --mem=4G --time=0-03:00 --array=${2:-1}-${3:-1000} \
      --export=ALL,MODE=power,MISS=${4:-0:0.25:0.5:0.75:0.9},RHO=${5:-0.75} "$0"
  fi
  exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
if [ "$MODE" = t1e ]; then
  Rscript simulation_code/rev1/true_geno_typeIerror_rho075_rev1.R $SLURM_ARRAY_TASK_ID $NARR ${MISS//:/,} "" ${OUTDIR:-}
else
  Rscript simulation_code/rev1/true_geno_power_rho075_rev1.R $SLURM_ARRAY_TASK_ID ${MISS//:/,} "" "${OUTDIR:-}" ${RHO:-0.75}
fi

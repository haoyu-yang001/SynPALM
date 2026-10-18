#!/bin/bash
# True-genotype type I error (all rho, calibrated autosomal SNPs; true_geno_typeIerror_calgeno_rev1.R):
#   ./submit_true_geno_t1e_calgeno_rev1.sh [first] [last] [n_arrays=40] [n_reps=17]
# Task t (1 .. n_reps * n_arrays) runs replicate (t-1)/n_arrays + 1 on SNP array (t-1)%n_arrays + 1:
# ~14,400 SNPs x 5 missing rates x 4 rho, ~45-60 min. Finished (replicate, array, miss)
# outputs are skipped, so requeued / resubmitted tasks resume.
#SBATCH --job-name=rev1_tg_t1e_cal
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/tgcal_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/tgcal_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  NA_=${3:-40}; NR=${4:-17}
  sbatch --array=${1:-1}-${2:-$((NA_ * NR))} --export=ALL,NARR=$NA_ "$0"; exit 0
fi
REP=$(( (SLURM_ARRAY_TASK_ID - 1) / NARR + 1 ))
ARR=$(( (SLURM_ARRAY_TASK_ID - 1) % NARR + 1 ))
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/true_geno_typeIerror_calgeno_rev1.R $REP $ARR $NARR

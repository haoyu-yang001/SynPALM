#!/bin/bash
# Variance components (variance_components_rev1.R) with each array split into blocks of
# replicates, for the settings whose full array exceeds the 3-day limit (missing rate 0.9):
#   ./submit_vc_split_rev1.sh <miss> [reps_per_task=2] [n_arrays=20] [n_rep=10]
# Task t runs SNP array (t-1)%n_arrays+1, replicates b*k+1..(b+1)*k with b=(t-1)/n_arrays;
# every (replicate, rho) fit is saved as it finishes, so preempted tasks resume, and the
# results equal those of an unsplit run (checked on a test array).
#SBATCH --job-name=rev1_vcsplit
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/vcsplit_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/vcsplit_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=2-00:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  K=${2:-2}; NA_=${3:-20}; NR=${4:-10}
  sbatch --job-name=rev1_vcsplit_$1 --array=1-$(( NA_ * NR / K )) --export=ALL,MISS=$1,K=$K,NARR=$NA_,NREP=$NR "$0"; exit 0
fi
A=$(( (SLURM_ARRAY_TASK_ID - 1) % NARR + 1 )); B=$(( (SLURM_ARRAY_TASK_ID - 1) / NARR ))
module load R/4.4.1-fasrc01
export OMP_NUM_THREADS=1
cd $BASE
VC_RHO_FIRST=${VC_RHO_FIRST-0.75} Rscript simulation_code/rev1/variance_components_rev1.R $MISS $A $NREP 10000 $(( B * K + 1 )) $(( (B + 1) * K ))

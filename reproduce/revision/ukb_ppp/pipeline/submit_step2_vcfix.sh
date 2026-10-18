#!/bin/bash
# Step 2 rerun after the 2026-10-13 VC fix, for proteins on lines
# <start>-<end> of ORDER_FILE whose step-1 cache exists.
#   ./submit_step2_vcfix.sh 1 <n> [max_parallel]
#SBATCH --job-name=cf3_s2fix
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_vcfix/s2fix_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_vcfix/s2fix_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=0-6:00
BASE=<PROJECT_DIR>
ORDER_FILE=${ORDER_FILE:-$BASE/code/sh_run/run_noPRS_noPC/order_vcfix_main.txt}
RES=$BASE/results/allukb_protein_pQTL_european_discovery
TAG=${SYNPALM_S1_TAG:-_cf_noPRS_noPC}
mkdir -p $BASE/logs/noPRS_noPC_vcfix
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=${1:-1}-${2:-2916}%${3:-400} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
# MKL in place of R's reference BLAS (sequential: the parallelism is the forked workers)
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
p=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ORDER_FILE" | tr -d '[:space:]')
name=$(head -1 data/df_allUKB_eur1.tsv | cut -f$((48 + p)))
[ -s "$RES/$name/${name}${TAG}_scan_allSNPs.rds" ] && { echo "skip $name"; rm -f $BASE/logs/noPRS_noPC_vcfix/s2fix_${SLURM_ARRAY_TASK_ID}.{out,err}; exit 0; }
[ -s "$RES/$name/${name}_synthetic_phenotype${TAG}.rds" ] || { echo "no step-1 cache for $name"; exit 3; }
export SYNPALM_S1_TAG=$TAG SYNPALM_VALIDATE=0
Rscript code/synpalm_main_code/step2_protein_fast.R $p 10 && rm -f $BASE/logs/noPRS_noPC_vcfix/s2fix_${SLURM_ARRAY_TASK_ID}.{out,err}

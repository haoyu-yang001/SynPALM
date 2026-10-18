#!/bin/bash
# Step 1 (cross-fitted synthetic phenotype, no PRS, no PCs) for proteins on
# lines <start>-<end> of ORDER_FILE.   ./submit_step1.sh 1 2916 [max_parallel]
#SBATCH --job-name=cf3_s1
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_step1/s1_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_step1/s1_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-3:00
BASE=<PROJECT_DIR>
ORDER_FILE=${ORDER_FILE:-$BASE/code/sh_run/order_cor_vec7_new.txt}
RES=$BASE/results/allukb_protein_pQTL_european_discovery
TAG=_cf_noPRS_noPC
mkdir -p $BASE/logs/noPRS_noPC_step1
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=${1:-1}-${2:-2916}%${3:-400} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
cd $BASE
export SYNPALM_PERMUTE_Y=${SYNPALM_PERMUTE_Y:-FALSE}
p=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ORDER_FILE" | tr -d '[:space:]')
name=$(head -1 data/df_allUKB_eur1.tsv | cut -f$((48 + p)))
suffix=$TAG; [ "$SYNPALM_PERMUTE_Y" = TRUE ] && suffix=${TAG}_PERMUTED${SYNPALM_PERM_SEED:-1}
if [ -s "$RES/$name/${name}_synthetic_phenotype${suffix}.rds" ]; then echo "skip $name"; rm -f $BASE/logs/noPRS_noPC_step1/s1_${SLURM_ARRAY_TASK_ID}.{out,err}; exit 0; fi
Rscript code/synpalm_main_code/step1_protein_crossfit.R $p && rm -f $BASE/logs/noPRS_noPC_step1/s1_${SLURM_ARRAY_TASK_ID}.{out,err}

#!/bin/bash
# Post-processing in parallel: task t of N handles lines with (line - 1) %% N == t - 1
# of the order file (skips proteins without a scan or already done). Then run
# aggregate_all.R once (submit_post.sh does both serially).
#   ./submit_post_array.sh [N=40]
#SBATCH --job-name=cf3_posta
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/posta_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/posta_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-12:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=1-${1:-40} --export=ALL,NTASK=${1:-40} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
cd $BASE
RES=results/allukb_protein_pQTL_european_discovery
H=$(head -1 data/df_allUKB_eur1.tsv)
for line in $(seq $SLURM_ARRAY_TASK_ID $NTASK 2916); do
  p=$(sed -n "${line}p" code/sh_run/order_cor_vec7_new.txt | tr -d '[:space:]')
  name=$(echo "$H" | cut -f$((48 + p)))
  [ -s "$RES/$name/${name}_cf_noPRS_noPC_scan_allSNPs.rds" ] || continue
  [ -s "results/table_results_cf_noPRS_noPC/${name}_loci.rds" ] && continue
  Rscript code/post_cf3/post_protein.R $line || echo "FAILED line $line ($name)"
done

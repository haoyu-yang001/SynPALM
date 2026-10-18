#!/bin/bash
# Post-processing (positions, counts, clumping) for every protein whose scan is
# finished and not yet post-processed; then the cross-protein summary.
#   sbatch submit_post.sh        (safe to rerun: done proteins are skipped)
#SBATCH --job-name=cf3_post
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/post_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/post_%j.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-12:00
BASE=<PROJECT_DIR>
module load R/4.4.1-fasrc01
cd $BASE
RES=results/allukb_protein_pQTL_european_discovery
for line in $(seq 1 2916); do
  p=$(sed -n "${line}p" code/sh_run/order_cor_vec7_new.txt | tr -d '[:space:]')
  name=$(head -1 data/df_allUKB_eur1.tsv | cut -f$((48 + p)))
  [ -s "$RES/$name/${name}_cf_noPRS_noPC_scan_allSNPs.rds" ] || continue
  [ -s "results/table_results_cf_noPRS_noPC/${name}_loci.rds" ] && continue
  Rscript code/post_cf3/post_protein.R $line || echo "FAILED line $line ($name)"
done
Rscript code/post_cf3/aggregate_all.R

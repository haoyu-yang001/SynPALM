#!/bin/bash
# Post-processing of all finished, not yet processed proteins with P parallel
# processes inside one job (for partitions with a per-user job limit, e.g. test).
#   sbatch -p test -c 8 --mem=64G -t 0-6:00 post_parallel_local.sh [P=8]
#SBATCH --job-name=cf3_postp
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/postp_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/postp_%j.err
BASE=<PROJECT_DIR>
module load R/4.4.1-fasrc01
cd $BASE
RES=results/allukb_protein_pQTL_european_discovery
H=$(head -1 data/df_allUKB_eur1.tsv)
todo=()
for line in $(seq 1 2916); do
  p=$(sed -n "${line}p" code/sh_run/order_cor_vec7_new.txt | tr -d '[:space:]')
  name=$(echo "$H" | cut -f$((48 + p)))
  [ -s "$RES/$name/${name}_cf_noPRS_noPC_scan_allSNPs.rds" ] || continue
  [ -s "results/table_results_cf_noPRS_noPC/${name}_loci.rds" ] && continue
  todo+=($line)
done
echo "proteins to process: ${#todo[@]}"
printf "%s\n" "${todo[@]}" | xargs -P ${1:-8} -I{} sh -c 'Rscript code/post_cf3/post_protein.R {} || echo "FAILED line {}"'

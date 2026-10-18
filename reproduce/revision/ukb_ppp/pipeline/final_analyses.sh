#!/bin/bash
# Cross-protein analyses of the no-PRS / no-PC rerun, in one job (test partition):
# aggregate_all.R, category_and_rho_strata.R, external_validation.R (P workers,
# then summarize) and the Fig 4 panels.
#   sbatch -p test -c 8 --mem=64G -t 0-4:00 final_analyses.sh [P=6]
#SBATCH --job-name=cf3_final
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/final_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/final_%j.err
BASE=<PROJECT_DIR>
module load R/4.4.1-fasrc01
cd $BASE
P=${1:-6}
echo "== aggregate_all";   Rscript code/post_cf3/aggregate_all.R || exit 1
echo "== rho strata";      Rscript code/post_cf3/category_and_rho_strata.R
echo "== external";        seq 1 $P | xargs -P $P -I{} Rscript code/post_cf3/external_validation.R {} $P
Rscript code/post_cf3/external_validation.R summarize
echo "== fig4";            Rscript code/post_cf3/figs/fig4_panels_cf.R
echo done

#!/bin/bash
# Final cross-protein analyses of the no-PRS / no-PC rerun, once all 2,916 scans
# exist (one job, test partition): post-processing of any unprocessed protein,
# aggregate_all.R, category_and_rho_strata.R, external_validation.R (P workers;
# per-protein outputs are kept, missing ones are added) + summarize,
# headline_numbers.R, SuSiE fine-mapping (missing proteins) + summarize, and the
# Fig 4 / Fig 5 panels.
#   sbatch -p test -c 32 --mem=192G -t 0-6:00 final_all.sh [P=32]
#SBATCH --job-name=cf3_all
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/all_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/all_%j.err
BASE=<PROJECT_DIR>
cd $BASE
P=${1:-32}
echo "== post";      bash code/sh_run/run_noPRS_noPC/post_parallel_local.sh 8
module load R/4.4.1-fasrc01
export OMP_NUM_THREADS=1
echo "== aggregate"; Rscript code/post_cf3/aggregate_all.R || exit 1
echo "== strata";    Rscript code/post_cf3/category_and_rho_strata.R
echo "== external";  seq 1 $P | xargs -P $P -I{} Rscript code/post_cf3/external_validation.R {} $P
Rscript code/post_cf3/external_validation.R summarize
echo "== headline";  Rscript code/post_cf3/headline_numbers.R
echo "== susie";     Rscript code/post_cf3/finemap_susie.R 1 1 $P && Rscript code/post_cf3/finemap_susie.R summarize
echo "== figs";      Rscript code/post_cf3/figs/fig4_panels_cf.R; Rscript code/post_cf3/figs/fig5b_circos_susie.R
echo done

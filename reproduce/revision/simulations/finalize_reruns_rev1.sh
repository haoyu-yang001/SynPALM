#!/bin/bash
# Final step after the reruns of the settings kept from the submission (power rho 0/0.25/0.5,
# bias rho 0, true-genotype power rho 0/0.25/0.5; 2026-10-16): check completeness and
# readability of every output, merge, redraw the figures, rebuild the composites and copy
# the figures used by the paper into paper/figs (local only; nothing is pushed).
#   sbatch -p test -c 4 --mem=32G -t 0-3:00 code/sh_run/finalize_reruns_rev1.sh   [FORCE=1 to skip the completeness stop]
#SBATCH -J finalize_rev1
#SBATCH -o <PROJECT_DIR>/logs/sim_rev1/finalize_%j.out
set -u
BASE=<PROJECT_DIR>
cd $BASE
module load R/4.4.1-fasrc01
R1=simulation_code/results_rev1; F=$R1/figs; P=paper/figs
step() { echo; echo "=== $(date +%H:%M:%S) $*"; }

step "1. completeness"
python3 code/sh_run/check_sim_outputs.py | grep -E "rho0|rho025|rho05|bias_rho" | tee /tmp/finalize_chk_$$.txt
if grep -qE "missing +[1-9]" /tmp/finalize_chk_$$.txt && [ "${FORCE:-0}" != 1 ]; then
  echo "reruns incomplete; stop (FORCE=1 to continue with the submitted results for the missing settings)"; exit 1
fi

step "2. every rerun file readable (truncated files from killed workers would break the merge)"
Rscript -e '
dirs <- file.path("simulation_code/results_rev1", c("power_rho000","power_rho025","power_rho050","bias_rho000",
                                                   "true_geno_power_rho000","true_geno_power_rho025","true_geno_power_rho050"))
bad <- unlist(lapply(dirs, function(d) { f <- list.files(d, full.names = TRUE)
  f[!vapply(f, function(x) tryCatch({ readRDS(x); TRUE }, error = function(e) FALSE), logical(1))] }))
cat("unreadable files:", length(bad), "\n"); if (length(bad)) { writeLines(bad); quit(status = 1) }' || exit 1

step "3. merge power / bias (replaces the links to the submitted rho <= 0.5 files)"
Rscript simulation_code/rev1/merge_power_bias_rev1.R || exit 1

step "4. figures"
Rscript simulation_code/rev1/plot_for_point_robustness_rev1.R || exit 1          # Fig 2A
Rscript simulation_code/rev1/power_analysis_kinsetting_plots_rev1.R || exit 1    # Fig 2C, Fig 3, S5-S7
Rscript simulation_code/rev1/analyze_true_geno_rev1.R || exit 1                  # Table S10, Figs S11-S12

step "5. composites (copied to paper/figs by the build script)"
bash code/sh_run/build_composite_figs.sh fig2 fig3

step "6. copy the single figures used by supp.tex"
for k in cousin halfsiblings siblings; do cp $F/power_Syn-PALM_ObsG_kinsetting_$k.pdf $P/power_Syn-PALM_ObsG_kinsetting_${k}_rev1.pdf; done
cp $F/power_Syn-PALM_ObsG_true_geno_rev1.pdf $F/power_ratio_Syn-PALM_vs_obsG_true_geno_rev1.pdf \
   $F/power_ratio_Syn-PALM_vs_SynSurr_true_geno_rev1.pdf $P/
ls -la --time-style=+%m-%d_%H:%M $P | grep -E "fig2_rev1|fig3_rev1|kinsetting|true_geno"

step "7. numbers quoted in the text (check against manuscript/supp/reply)"
Rscript -e '
x <- read.csv("simulation_code/results_rev1/figs/power_true_geno_rev1.csv"); x <- x[abs(x$pve_g - 5e-4) < 1e-12, ]
x$ratio <- x$power_SynSurrG / x$power_ObsG; x$ratio_ss <- x$power_SynSurrG / x$power_SynSurr
cat("true-genotype ratio vs Obs-GWAS (h2 0.05%), rho = 0 or miss = 0:", sprintf("%.3f", range(x$ratio[x$rho == 0 | x$miss == 0])), "\n")
cat("max ratio (miss 0.9, rho 0.75):", sprintf("%.2f", x$ratio[x$miss == 0.9 & x$rho == 0.75]), "; vs SynSurr range", sprintf("%.2f", range(x$ratio_ss)), "\n")
print(x[order(x$rho, x$miss), c("miss", "rho", "power_SynSurrG", "power_ObsG", "ratio")], row.names = FALSE)'
echo; echo "done $(date)"

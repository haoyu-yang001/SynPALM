## Reply to Reviewer 1, Comment 7 and minor comment (e) (revision, 2026-10-19): number of
## primary associations of Syn-PALM and Obs-GWAS
##   (i)  at the proteome-wide Bonferroni threshold 5e-8 / 2,916 = 1.71e-11, clumped as in
##        post_protein.R (+-1 Mb windows around significant SNPs merged, sentinel = lowest P,
##        MHC chr6:25.5-34.0 Mb removed); the thresholds 5e-8 ... 1e-10 come from
##        aggregate_all.R (summary_per_protein.csv), and are recomputed here as a check;
##   (ii) at 5e-8 with +-250 kb and +-500 kb windows instead of +-1 Mb.
## Reads the per-protein web-portal summary statistics
## (results/web_SynGWAS_summary_results_cf_noPRS_noPC/).
##   sbatch -p test -c 16 --mem=32G -t 0-2:00 --wrap \
##     "module load R/4.4.1-fasrc01; Rscript code/post_cf3/threshold_window_sensitivity.R"
## Writes results/table_results_cf_noPRS_noPC/threshold_window_sensitivity{_per_protein,}.csv
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(parallel) })
setwd("<PROJECT_DIR>")
WEB <- "results/web_SynGWAS_summary_results_cf_noPRS_noPC/"
OUT <- "results/table_results_cf_noPRS_noPC/"
NC  <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
files <- list.files(WEB, "_negative_log10_pvals_v2geno_crossfit\\.rds$", full.names = TRUE)
grid <- rbind(data.table(thr = c(5e-8, 1e-8, 1e-9, 1e-10, 5e-8 / 2916), win = 1e6),
              data.table(thr = 5e-8, win = c(2.5e5, 5e5)))
n_loci <- function(s, lp, thr, win) {
  s <- s[s[[lp]] >= -log10(thr)]
  if (!nrow(s)) return(0L)
  length(reduce(GRanges(as.character(s$chr), IRanges(s$POS38 - win, s$POS38 + win))))
}
res <- rbindlist(mclapply(files, function(f) {
  x <- readRDS(f)
  x <- x[!(chr == 6 & POS38 >= 25.5e6 & POS38 <= 34e6)]
  m <- -log10(5e-8)
  x <- x[SynSurrG_negative_log10_pvalues >= m | ObsG_negative_log10_pvalues >= m]
  p <- sub("_negative_log10_pvals_v2geno_crossfit\\.rds$", "", basename(f))
  grid[, .(protein = p, thr, win,
           SynG = mapply(function(t, w) n_loci(x, "SynSurrG_negative_log10_pvalues", t, w), thr, win),
           ObsG = mapply(function(t, w) n_loci(x, "ObsG_negative_log10_pvalues", t, w), thr, win))]
}, mc.cores = NC))
fwrite(res, paste0(OUT, "threshold_window_sensitivity_per_protein.csv"))
tot <- res[, .(proteins = uniqueN(protein), SynG = sum(SynG), ObsG = sum(ObsG),
               proteins_with_SynG = sum(SynG > 0), proteins_with_ObsG = sum(ObsG > 0)), by = .(thr, win)]
tot[, ratio := SynG / ObsG]
## check against aggregate_all.R at +-1 Mb
per <- fread(paste0(OUT, "summary_per_protein.csv"))
for (t in c("5e-08", "1e-08", "1e-09", "1e-10"))
  cat(sprintf("check %s: aggregate SynG %d ObsG %d\n", t, sum(per[[paste0("SynG_loci_", t)]]), sum(per[[paste0("ObsG_loci_", t)]])))
fwrite(tot, paste0(OUT, "threshold_window_sensitivity.csv"))
print(tot)

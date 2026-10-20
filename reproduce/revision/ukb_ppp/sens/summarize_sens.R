## Summary of the sensitivity analysis (revision 2026-10-19; reply to Reviewer 1,
## Comments 4 and 8): for each protein and setting, compared with the main analysis on
## the same SNP set (main = results/web_SynGWAS_summary_results_cf_noPRS_noPC/):
##   rho_oof (main, new); correlation of the Syn-PALM -log10 P over the 20,000 random SNPs;
##   genome-wide significant SNPs (P < 5e-8) and primary associations (+-1 Mb merged
##   windows, MHC removed) of the main and the new run and their overlap; the change of
##   the Syn-PALM estimate at the main GWS SNPs in units of its main SE; and, as a check,
##   the largest difference of the Obs-GWAS -log10 P (Obs-GWAS does not use Yhat).
##   Rscript code/post_cf3/sens/summarize_sens.R
## Writes results/sens_cf_noPRS_noPC/sens_per_protein.csv and sens_summary.csv
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
setwd("<PROJECT_DIR>")
OUT <- "results/sens_cf_noPRS_noPC/"
pr  <- fread(paste0(OUT, "proteins.csv"))
rnd <- readLines(paste0(OUT, "snps_random20k.txt"))
tags <- c(fold_seed_1 = "_fs20261101", fold_seed_2 = "_fs20261102", predictors_50 = "_np50", predictors_200 = "_np200")
THR <- -log10(5e-8)
loci <- function(d, lp) {
  s <- d[get(lp) >= THR & !(chr == 6 & POS38 >= 25.5e6 & POS38 <= 34e6)]
  if (!nrow(s)) return(GRanges())
  reduce(GRanges(as.character(s$chr), IRanges(s$POS38 - 1e6, s$POS38 + 1e6)))
}
rows <- list()
for (i in seq_len(nrow(pr))) {
  p <- pr$protein[i]; D <- sprintf("results/allukb_protein_pQTL_european_discovery/%s/", p)
  w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
  for (s in names(tags)) {
    f <- sprintf("%s%s_cf_noPRS_noPC%s_scan_subset.rds", D, p, tags[[s]])
    if (!file.exists(f)) next
    x <- readRDS(f); s1 <- readRDS(sprintf("%s%s_synthetic_phenotype_cf_noPRS_noPC%s.rds", D, p, tags[[s]]))
    d <- merge(w[, .(rsid, chr, POS38, m_lp = SynSurrG_negative_log10_pvalues, m_b = SynSurrG_hat_beta,
                     m_se = sqrt(SynSurrG_var_hat_beta), mo_lp = ObsG_negative_log10_pvalues)],
               data.table(rsid = rownames(x), n_lp = x[, "SynSurrG_negative_log10_pvalues"],
                          n_b = x[, "SynSurrG_hat_beta"], no_lp = x[, "ObsG_negative_log10_pvalues"]), by = "rsid")
    gm <- d[m_lp >= THR, rsid]; gn <- d[n_lp >= THR, rsid]
    lm <- loci(d, "m_lp"); ln <- loci(d, "n_lp")
    rows[[length(rows) + 1]] <- data.table(
      protein = p, q = pr$q[i], setting = s, rho_main = pr$rho_oof[i], rho_new = s1$rho_oof,
      r_nlp_random = d[rsid %in% rnd, cor(m_lp, n_lp)],
      gws_main = length(gm), gws_new = length(gn), gws_both = length(intersect(gm, gn)),
      loci_main = length(lm), loci_new = length(ln),
      loci_main_found = if (length(lm)) sum(overlapsAny(lm, ln)) else 0L,
      loci_new_found = if (length(ln)) sum(overlapsAny(ln, lm)) else 0L,
      dz_med = d[m_lp >= THR, median(abs(n_b - m_b) / m_se)], dz_max = d[m_lp >= THR, max(abs(n_b - m_b) / m_se)],
      obs_maxdiff = d[, max(abs(no_lp - mo_lp))])
  }
}
r <- rbindlist(rows)
fwrite(r, paste0(OUT, "sens_per_protein.csv"))
r[, type := fifelse(grepl("fold", setting), "fold seed", setting)]
sm <- r[, .(runs = .N, rho_diff_median = median(rho_new - rho_main), rho_diff_range = sprintf("%.3f to %.3f", min(rho_new - rho_main), max(rho_new - rho_main)),
            r_nlp_median = median(r_nlp_random), r_nlp_min = min(r_nlp_random),
            gws_main = sum(gws_main), gws_new = sum(gws_new), gws_overlap = sum(gws_both) / sum(gws_main),
            loci_main = sum(loci_main), loci_new = sum(loci_new),
            loci_main_recovered = sum(loci_main_found) / sum(loci_main), loci_new_in_main = sum(loci_new_found) / sum(loci_new),
            dz_median = median(dz_med), dz_max = max(dz_max), obs_check = max(obs_maxdiff)), by = type]
fwrite(sm, paste0(OUT, "sens_summary.csv"))
print(r, digits = 3); print(sm, digits = 3)

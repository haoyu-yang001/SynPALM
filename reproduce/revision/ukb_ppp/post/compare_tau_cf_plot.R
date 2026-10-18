## Supplementary figure "twoSGWAS" (revision): merges the per-array output of
## code/post_cf3/compare_tau_cf.R (refit mode) and plots, for the 2,000 SNPs most
## associated with the synthetic phenotype, Syn-PALM -log10 P with variance
## components from the global null (alpha_G = 0, production) against -log10 P
## with SNP-specific variance components (alpha_G free), with the identity line,
## in the style of the submitted figure (code/combine_plot_compare_taus.R).
## Also checks the global-null column against the production genome-wide scan
## (<p>_cf_noPRS_noPC_scan_allSNPs.rds) and writes agreement statistics.
##   Rscript code/post_cf3/compare_tau_cf_plot.R [n_array=100]
## env: SYNPALM_CT_PROTEIN protein index (default 1571 = lep)
## Writes:
##   results/compare_tau_cf_noPRS_noPC/<p>_TOTAL_compare_tau_cf_noPRS_noPC.rds  merged table
##     (per-array files in temp/ are removed once the merged file is saved)
##   results/compare_tau_cf_noPRS_noPC/<p>_compare_tau_summary.csv               agreement statistics
##   results/figs_cf_noPRS_noPC/SynGWASalphaG_rev1.pdf                           figure
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
NARR <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 100L
protein_index <- as.integer(Sys.getenv("SYNPALM_CT_PROTEIN", "1571"))
protein_names <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
p   <- protein_names[protein_index]
OUT <- "results/compare_tau_cf_noPRS_noPC/"
FIG <- "results/figs_cf_noPRS_noPC/"; dir.create(FIG, showWarnings = FALSE)
tot_file <- paste0(OUT, p, "_TOTAL_compare_tau_cf_noPRS_noPC.rds")

if (!file.exists(tot_file)) {
  fs <- sprintf("%stemp/%s_refit_arr%dof%d.rds", OUT, p, seq_len(NARR), NARR)
  miss <- which(!file.exists(fs))
  if (length(miss)) stop(length(miss), " array output(s) missing: ", paste(head(miss, 20), collapse = ","))
  L   <- lapply(fs, readRDS)
  vc0 <- attr(L[[1]], "vc_global")
  stopifnot(all(vapply(L, function(x) isTRUE(all.equal(attr(x, "vc_global"), vc0)), logical(1))))
  res <- do.call(rbind, L)
  top <- readRDS(paste0(OUT, p, "_top2000_yhat_snps.rds"))
  stopifnot(identical(rownames(res), top$rsid))
  res <- data.table(top[, .(rsid, chr, POS38, yhat_beta, yhat_se, yhat_log10p)], as.data.table(res))
  attr(res, "vc_global") <- vc0
  saveRDS(res, tot_file)
  invisible(file.remove(fs))
}
res <- readRDS(tot_file)
vc0 <- attr(res, "vc_global")

## global-null column vs the production genome-wide scan
prod <- readRDS(sprintf("results/allukb_protein_pQTL_european_discovery/%s/%s_cf_noPRS_noPC_scan_allSNPs.rds", p, p))
d_prod <- max(abs(prod[res$rsid, "SynSurrG_negative_log10_pvalues"] - res$SynSurrG_negative_log10_pvalues))
cat(sprintf("global-null -log10 P vs production scan: max |diff| = %.3g\n", d_prod))

x <- res$SynSurrG_negative_log10_pvalues; y <- res$SynSurrG_each_negative_log10_pvalues
thr <- -log10(5e-8); ad <- abs(x - y)
S <- data.table(
  protein = p, n_snps = nrow(res), n_chr = uniqueN(res$chr),
  yhat_log10p_min = min(res$yhat_log10p), yhat_log10p_max = max(res$yhat_log10p),
  pearson = cor(x, y), spearman = cor(x, y, method = "spearman"),
  max_abs_diff = max(ad), median_abs_diff = median(ad), mean_abs_diff = mean(ad),
  max_rel_diff = max(ad / pmax(x, 1)), rsid_max_diff = res$rsid[which.max(ad)],
  slope_through_origin = sum(x * y) / sum(x^2),
  n_sig_global = sum(x > thr), n_sig_each = sum(y > thr),
  n_sig_global_only = sum(x > thr & y <= thr), n_sig_each_only = sum(y > thr & x <= thr),
  max_log10p_global = max(x), max_log10p_each = max(y),
  tau_S2_global = vc0[["tau_S2"]], tau_S2_each_min = min(res$each_tau_S2), tau_S2_each_max = max(res$each_tau_S2),
  tau_TS_global = vc0[["tau_TS"]], tau_TS_each_min = min(res$each_tau_TS), tau_TS_each_max = max(res$each_tau_TS),
  n_projected_G_each = sum(res$each_projected_G), projected_G_global = vc0[["projected_G"]],
  median_secs_per_refit = median(res$secs), maxdiff_vs_production = d_prod)
fwrite(S, paste0(OUT, p, "_compare_tau_summary.csv"))
print(t(S))

lim <- c(0, max(x, y) * 1.03)
p1 <- ggplot(res, aes(x = SynSurrG_negative_log10_pvalues, y = SynSurrG_each_negative_log10_pvalues)) +
  geom_point(alpha = 1, color = "#4575b4", size = 1) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color = "#d13027", linetype = "dashed") +
  labs(title = "Comparison of -log10(P-values)",
       x = bquote("Syn-PALM (" * alpha[G] == 0 * ") -log"[10] * "(p)"),
       y = bquote("Syn-PALM (" * alpha[G] != 0 * ") -log"[10] * "(p)")) +
  scale_x_continuous(limits = lim) + scale_y_continuous(limits = lim) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  coord_fixed()
ggsave(plot = p1, filename = paste0(FIG, "SynGWASalphaG_rev1.pdf"),
       device = "pdf", height = 3.7, width = 3.7, units = "in")
cat("written:", paste0(FIG, "SynGWASalphaG_rev1.pdf"), "\n")

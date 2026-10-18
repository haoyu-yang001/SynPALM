## Fig S16 (revision): QQ plots of Syn-PALM (and Obs-GWAS) under the null for
## the 10 proteins of the submitted figure, after the joint permutation of
## (Y, Yhat) done by step2_protein_fast_jointperm.R.
##   Rscript code/post_cf3/qq_permuted.R [seed=1]
## Writes results/figs_cf_noPRS_noPC/qq_ukbPPP_permuted.{png,pdf} and a lambda_GC table.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
setwd("<PROJECT_DIR>")
SEED <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "1"
PROT <- c("ccl24", "chit1", "cst3", "enpp5", "il18r1", "lilrb2", "lrpap1", "msmb", "pilra", "tnfrsf10c")
OUT  <- "results/figs_cf_noPRS_noPC/"; dir.create(OUT, showWarnings = FALSE)
R    <- "results/allukb_protein_pQTL_european_discovery/"
ann  <- readRDS("data/snp_annot_v2geno.rds")              # autosomes only, as in the paper

thin <- function(lp) {             # keep every SNP in the tail, thin the bulk
  n <- length(lp); o <- sort(lp, decreasing = TRUE)
  e <- -log10((seq_len(n) - 0.5) / n)
  i <- unique(c(seq_len(min(n, 5000)), round(seq(1, n, length.out = 20000))))
  lo <- -log10(qbeta(0.975, i, n - i + 1)); hi <- -log10(qbeta(0.025, i, n - i + 1))
  data.table(expected = e[i], observed = o[i], lo = lo, hi = hi)
}
pts <- list(); lam <- list()
for (p in PROT) {
  f <- sprintf("%s%s/%s_cf_noPRS_noPC_jointPERM%s_scan_allSNPs.rds", R, p, p, SEED)
  if (!file.exists(f)) { cat("missing:", p, "\n"); next }
  x <- readRDS(f)
  x <- x[rownames(x) %in% ann$rsid[ann$chr %in% as.character(1:22)], , drop = FALSE]
  for (m in c("SynSurrG", "ObsG")) {
    lp <- x[, paste0(m, "_negative_log10_pvalues")]; lp <- lp[is.finite(lp)]
    chi <- qchisq(10^(-lp), 1, lower.tail = FALSE)
    lam[[length(lam) + 1]] <- data.table(protein = p, method = m, n_snps = length(lp),
                                         lambda_gc = median(chi) / qchisq(0.5, 1),
                                         n_p_lt_5e8 = sum(lp >= -log10(5e-8)),
                                         max_log10p = max(lp))
    pts[[length(pts) + 1]] <- cbind(protein = p, method = m, thin(lp))
  }
}
lam <- rbindlist(lam); pts <- rbindlist(pts)
fwrite(lam, paste0(OUT, "qq_ukbPPP_permuted_lambda.csv"))
print(lam)
pts[, method := factor(method, c("SynSurrG", "ObsG"), c("Syn-PALM", "Obs-GWAS"))]
lab <- lam[method == "SynSurrG", .(protein, lab = sprintf("lambda == %.3f", lambda_gc))]
mx <- max(pts$observed, pts$hi) * 1.05
g <- ggplot(pts[method == "Syn-PALM"], aes(expected, observed)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  geom_point(size = 0.5, colour = "#0072B5") +
  geom_text(data = lab, aes(x = 0.2, y = mx * 0.92, label = lab), parse = TRUE,
            hjust = 0, size = 3) +
  facet_wrap(~ toupper(protein), nrow = 2) +
  coord_cartesian(ylim = c(0, mx)) +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P)))) +
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank())
ggsave(paste0(OUT, "qq_ukbPPP_permuted.png"), g, width = 11, height = 5, dpi = 300)
ggsave(paste0(OUT, "qq_ukbPPP_permuted.pdf"), g, width = 11, height = 5)
cat("wrote", paste0(OUT, "qq_ukbPPP_permuted.png"), "\n")

## =====================================================================
## marginal_truth_semisyn.R -- in-sample marginal effect of every design SNP
## on the planted part of Y* (and on the planted part of each surrogate).
##
## The planted causal SNPs are >= 1-2 Mb from every other planted SNP, but a
## SNP can still be in weak long-range LD with a causal SNP, in which case its
## marginal effect on Y* is not exactly 0. Because the planted scores are fixed
## across replicates, the marginal truth of SNP j is
##     t_j = cov_n(G_j, score_Y) / var_n(G_j)          (analysis cohort, n = 402,086)
## and similarly a_j(col) for the surrogate columns. analyze_semisyn.R uses t_j
## as the truth for coverage/bias and flags null SNPs with a non-negligible t_j.
## Also stores the largest r^2 of each SNP with a causal SNP on its chromosome.
## Output: simulation_code/results_rev1/semisyn/design_truth.rds (per-SNP, summary level)
##   Rscript marginal_truth_semisyn.R
## =====================================================================
setwd("<PROJECT_DIR>")
suppressMessages({library(data.table)})
OUT <- "simulation_code/results_rev1/semisyn/"
design <- readRDS(paste0(OUT, "design.rds")); snps <- design$snps
DI <- readRDS(paste0(OUT, "design_indiv.rds"))
FD <- readRDS("data/crossfit_folds_K5.rds")
G0 <- BEDMatrix::BEDMatrix("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
rows <- match(unname(as.character(FD$analysis_eid)), rownames(G0)); stopifnot(!anyNA(rows))
sc <- scale(DI$score, scale = FALSE)
causal <- snps$rsid[snps$cat %in% c("causal_Y", "causal_YS")]
Gc <- as.matrix(G0[rows, causal]); Gc[is.na(Gc)] <- 0; Gc <- scale(Gc)
t0 <- Sys.time()
res <- parallel::mclapply(split(snps$rsid, ceiling(seq_along(snps$rsid) / 200)), function(ids) {
  G <- as.matrix(G0[rows, ids, drop = FALSE]); G[is.na(G)] <- 0
  G <- sweep(G, 2, colMeans(G)); v <- colSums(G^2)
  eff <- crossprod(G, sc) / v                                    # m x 6 marginal slopes
  r <- crossprod(G / sqrt(v / (nrow(G) - 1)), Gc) / (nrow(G) - 1) # correlations with causal SNPs
  chr_j <- snps$chr[match(ids, snps$rsid)]; chr_c <- snps$chr[match(causal, snps$rsid)]
  r[outer(chr_j, chr_c, "!=")] <- 0
  self <- outer(ids, causal, "=="); r[self] <- 0
  data.frame(rsid = ids, t_Y = eff[, "Y"], eff[, colnames(eff) != "Y", drop = FALSE],
             max_r2_causal = apply(r^2, 1, max), check.names = FALSE)
}, mc.cores = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))), mc.preschedule = FALSE)
stopifnot(all(vapply(res, is.data.frame, logical(1))))
tr <- do.call(rbind, res); rownames(tr) <- NULL
stopifnot(identical(tr$rsid, snps$rsid))
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
print(aggregate(cbind(abs_tY = abs(t_Y), max_r2_causal) ~ cat, data = cbind(tr, cat = snps$cat), FUN = function(x) c(mean = mean(x), max = max(x))))
saveRDS(tr, paste0(OUT, "design_truth.rds"))

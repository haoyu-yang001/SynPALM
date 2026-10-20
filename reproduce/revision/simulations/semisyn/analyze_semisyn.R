## =====================================================================
## analyze_semisyn.R -- summary of the semi-synthetic simulation
## (run_semisyn_rep.R outputs in simulation_code/results_rev1/semisyn/reps/).
##
## Truth on the analysis scale: the pipeline inverse-normal transforms Y_obs,
## so the planted effect b_j (per counted allele, Var(Y*) ~ 1) becomes
## kappa_r * b_j, kappa_r = mean d INT(y)/dy over the labelled Y* of replicate r
## (stored by run_semisyn_rep.R; ~ 1/sd(Y*) for a near-normal Y*).
##
## Reports, for Syn-PALM (SynSurrG), Obs-GWAS (ObsG) and the naive OLS GWAS of
## Yhat (contrast):
##   * type I error at 1e-3, 1e-4, 5e-8 for random nulls and S-acting nulls
##     (by column type and effect size), lambda_GC for random nulls;
##   * coverage of the 95% Wald CI, mean bias and slope of beta_hat on the
##     truth for the planted causal SNPs;
##   * power at 5e-8 and 1e-4 and power ratio Syn-PALM / Obs-GWAS
##     (replicate bootstrap 95% CI), variance ratio (efficiency);
##   * across-replicate variance of beta_hat (per SNP, pooled) vs the mean
##     reported variance, and MSE around the truth vs mean reported variance.
##
##   Rscript analyze_semisyn.R [subdir=reps]
## Output: simulation_code/results_rev1/semisyn/summary_semisyn_<subdir>.{csv,rds}
## =====================================================================
setwd("<PROJECT_DIR>")
suppressMessages(library(data.table))
OUT  <- "simulation_code/results_rev1/semisyn/"
sub  <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "reps"
fs   <- sort(list.files(paste0(OUT, sub), "^rep_\\d+\\.rds$", full.names = TRUE))
design <- readRDS(paste0(OUT, "design.rds")); snps <- design$snps
R <- length(fs); cat("replicates:", R, "\n")
reps <- lapply(fs, readRDS)
stopifnot(all(vapply(reps, function(x) identical(rownames(x), snps$rsid), logical(1))))
M <- c(SynPALM = "SynPALM", ObsG = "ObsG", YhatOLS = "YhatOLS")
arr <- function(col) sapply(reps, function(x) x[, col])                 # nSNP x R
P  <- lapply(M, function(m) arr(paste0(m, "_mlog10p")))
B  <- lapply(M, function(m) arr(paste0(m, "_beta")))
V  <- lapply(M, function(m) arr(paste0(m, "_var")))
kappa <- vapply(reps, attr, 0, "kappa")
rho   <- t(sapply(reps, function(x) attr(x, "rho")))
## truth = planted b_j (per counted allele) on the INT scale. The in-sample
## association of a SNP with the FIXED parts of Y* (h(S0), planted scores) is a
## per-SNP offset that is treated as sampling variation, as it would be in a
## new cohort. marginal_truth_semisyn.R gives r^2 of each SNP with the causal
## SNPs on its chromosome: S_null SNPs with r^2 >= 0.01 carry real Y signal
## (long-range LD) and are excluded in the "_clean" sets.
truth <- outer(snps$b_Y, kappa)                                          # nSNP x R
TF <- paste0(OUT, "design_truth.rds")
DT <- readRDS(TF); stopifnot(identical(DT$rsid, snps$rsid))
leak <- DT$max_r2_causal >= 0.01
cat("null-type SNPs with r2 >= 0.01 to a causal SNP:", sum(leak & snps$cat %in% c("null", "S_null")), "\n")
print(cbind(snps[leak & snps$cat == "S_null", c("rsid", "chr", "pos", "col")], r2 = DT$max_r2_causal[leak & snps$cat == "S_null"]))
info <- data.frame(rep = sapply(reps, function(x) attr(x, "settings")$rep),
                   seed = sapply(reps, function(x) attr(x, "settings")$seed),
                   rho_oof = rho[, "oof"], kappa = kappa,
                   time_min = sapply(reps, attr, "time_sec") / 60)
cat("rho_oof: mean", round(mean(info$rho_oof), 3), "range", paste(round(range(info$rho_oof), 3), collapse = "-"),
    "| kappa mean", round(mean(kappa), 4), "| minutes per rep", round(mean(info$time_min), 1), "\n")
rows <- list(); set.seed(1)
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(...)
sets <- list(null = snps$cat == "null",
             S_null = snps$cat == "S_null",
             S_null_hcol = snps$cat == "S_null" & snps$col != "plt",
             S_null_Z = snps$cat == "S_null" & snps$col == "plt",
             S_null_large = snps$cat == "S_null" & snps$qS >= 5e-3,
             S_null_clean = snps$cat == "S_null" & !leak)
for (s in names(sets)) for (m in names(M)) {
  p <- P[[m]][sets[[s]], , drop = FALSE]; N <- length(p)
  for (a in c(1e-3, 1e-4, 5e-8)) {
    ## SNPs are fixed across replicates, so the 95% CI resamples SNPs (rows)
    kj <- rowSums(p > -log10(a)); k <- sum(kj)
    bt <- replicate(2000, sum(kj[sample.int(length(kj), replace = TRUE)])) / N
    add(metric = "typeI", set = s, method = m, level = a, N = N, value = k / N,
        mcse = sqrt(a * (1 - a) / N), ratio_to_nominal = k / N / a,
        ci_lo = quantile(bt, 0.025) / a, ci_hi = quantile(bt, 0.975) / a,
        snps_hit = sum(kj > 0), max_hits_per_snp = max(kj))
  }
  chi <- qchisq(10^(-p), 1, lower.tail = FALSE)
  add(metric = "lambdaGC", set = s, method = m, level = NA, N = N, value = median(chi) / qchisq(0.5, 1),
      mcse = NA, ratio_to_nominal = NA)
}
csets <- list(causal_all = snps$cat %in% c("causal_Y", "causal_YS"),
              causal_Y = snps$cat == "causal_Y", causal_YS = snps$cat == "causal_YS")
for (q in sort(unique(snps$q[snps$q > 0]))) csets[[sprintf("causal_q%.0e", q)]] <- snps$q == q
boot_ratio <- function(num, den, B0 = 2000) {
  ## replicate bootstrap of sum(num)/sum(den), num/den: nSNP x R logical/numeric
  rn <- colSums(num); rd <- colSums(den)
  bt <- replicate(B0, { i <- sample.int(length(rn), replace = TRUE); sum(rn[i]) / sum(rd[i]) })
  quantile(bt, c(0.025, 0.975), na.rm = TRUE)
}
for (s in names(csets)) {
  ii <- csets[[s]]
  for (m in c("SynPALM", "ObsG")) {
    b <- B[[m]][ii, , drop = FALSE]; v <- V[[m]][ii, , drop = FALSE]; tr <- truth[ii, , drop = FALSE]
    cover <- abs(b - tr) <= qnorm(0.975) * sqrt(v)
    add(metric = "coverage95", set = s, method = m, level = NA, N = length(b), value = mean(cover),
        mcse = sqrt(0.95 * 0.05 / length(b)), ratio_to_nominal = mean(cover) / 0.95)
    add(metric = "bias_mean_over_meanSE", set = s, method = m, level = NA, N = length(b),
        value = mean(sign(tr) * (b - tr)) / mean(sqrt(v)), mcse = NA, ratio_to_nominal = NA)
    add(metric = "slope_beta_on_truth", set = s, method = m, level = NA, N = length(b),
        value = sum(b * tr) / sum(tr^2), mcse = NA, ratio_to_nominal = NA)
    for (a in c(5e-8, 1e-4)) add(metric = "power", set = s, method = m, level = a, N = length(b),
                                 value = mean(P[[m]][ii, ] > -log10(a)), mcse = NA, ratio_to_nominal = NA)
  }
  for (a in c(5e-8, 1e-4)) {
    hs <- P$SynPALM[ii, , drop = FALSE] > -log10(a); ho <- P$ObsG[ii, , drop = FALSE] > -log10(a)
    ci <- boot_ratio(hs, ho)
    add(metric = "power_ratio_Syn_over_Obs", set = s, method = "SynPALM/ObsG", level = a, N = length(hs),
        value = sum(hs) / sum(ho), mcse = NA, ratio_to_nominal = NA, ci_lo = ci[1], ci_hi = ci[2])
  }
  add(metric = "efficiency_meanVarObs_over_meanVarSyn", set = s, method = "SynPALM/ObsG", level = NA,
      N = sum(ii) * R, value = mean(V$ObsG[ii, ]) / mean(V$SynPALM[ii, ]), mcse = NA, ratio_to_nominal = NA)
}
## across-replicate variability vs reported SE (fixed truth for every SNP)
vsets <- c(list(causal = csets$causal_all), sets[c("S_null", "null")])
for (s in names(vsets)) for (m in c("SynPALM", "ObsG")) {
  ii <- vsets[[s]]; b <- B[[m]][ii, , drop = FALSE]; v <- V[[m]][ii, , drop = FALSE]; tr <- truth[ii, , drop = FALSE]
  if (R >= 2) add(metric = "SD_across_reps_over_meanSE", set = s, method = m, level = NA, N = length(b),
                  value = sqrt(mean(apply(b - tr, 1, var)) / mean(v)), mcse = NA, ratio_to_nominal = NA)
  add(metric = "RMSE_over_meanSE", set = s, method = m, level = NA, N = length(b),
      value = sqrt(mean((b - tr)^2) / mean(v)), mcse = NA, ratio_to_nominal = NA)
  add(metric = "mean_z2", set = s, method = m, level = NA, N = length(b),
      value = mean((b - tr)^2 / v), mcse = NA, ratio_to_nominal = NA)
}
tab <- do.call(rbind, lapply(rows, function(d) { for (c in c("ci_lo", "ci_hi", "snps_hit", "max_hits_per_snp")) if (is.null(d[[c]])) d[[c]] <- NA; d }))
rownames(tab) <- NULL
fwrite(tab, sprintf("%ssummary_semisyn_%s.csv", OUT, sub))
saveRDS(list(table = tab, info = info, R = R, leak = snps$rsid[leak]), sprintf("%ssummary_semisyn_%s.rds", OUT, sub))
fwrite(info, sprintf("%sreplicate_info_%s.csv", OUT, sub))
options(width = 200)
print(tab[, c("metric", "set", "method", "level", "N", "value", "ratio_to_nominal", "ci_lo", "ci_hi", "snps_hit", "max_hits_per_snp")], digits = 4)

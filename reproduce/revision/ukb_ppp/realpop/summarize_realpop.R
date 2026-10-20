## Summary of the real-data POP-GWAS and labelled-only analyses (revision 2026-10-20;
## reply to Reviewer 1, Comments 1 and 3, Reviewer 3 minor 1), from the outputs of
## step2_realpop.R and the main scans of the same 24 proteins.
##  (A) POP-GWAS (Miao et al.; pop_gwas/Fun.R, cal_qt_opt_vec, c = r N / (N + n)) on the
##      UKB-PPP data: labelled = individuals with the measured protein, unlabelled = the
##      rest; Z1/Z2 = GWAS of Yhat in the labelled / unlabelled sample, Z3 = GWAS of Y in
##      the labelled sample; r = corr(Y, Yhat | covariates) among the labelled.
##        POP_mm  mixed-model GWAS with the GRM restricted to each sample and effective
##                sample sizes n = n_lab var_OLS / var_MM (as in the simulations,
##                simulation_code/rev1/typeIerror_split_rev1.R)
##        POP_lm  ordinary least squares GWAS, n = n_lab, N = n_unlab
##      lambda_GC on the 20,000 random SNPs, genome-wide significant SNPs (P < 5e-8) and
##      primary associations (+-1 Mb merged windows, MHC removed) on the protein's SNP
##      set, for the real data and after the joint permutation of (Y, Yhat) (null).
##      SynSurr (unrelated subset) and Obs (unrelated subset) counts are genome-wide,
##      from the main scans (<p>_cf_noPRS_noPC_scan_allSNPs.rds, restricted to the
##      autosomal SNPs of the main summary results).
##  (B) Labelled-only Syn-PALM: SE ratio Syn-PALM / Obs-GWAS at the random SNPs (median),
##      GWS SNPs and primary associations, against the full analysis.
##   Rscript code/post_cf3/realpop/summarize_realpop.R
## Writes results/realpop_cf_noPRS_noPC/{pop_per_protein,pop_summary,labonly_per_protein,
##   labonly_summary,synsurr_unrelated_per_protein}.csv
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
setwd("<PROJECT_DIR>")
source("pop_gwas/Fun.R")
OUT <- "results/realpop_cf_noPRS_noPC/"
pr  <- fread("results/sens_cf_noPRS_noPC/proteins.csv")
rnd <- readLines("results/sens_cf_noPRS_noPC/snps_random20k.txt")
THR <- -log10(5e-8)
lp_to_chi <- function(lp) qchisq(-lp * log(10), 1, lower.tail = FALSE, log.p = TRUE)
lam <- function(lp) median(lp_to_chi(lp[is.finite(lp)])) / qchisq(0.5, 1)
z_to_lp <- function(z) -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
loci <- function(d, lp) {
  s <- d[d[[lp]] >= THR & !(chr == 6 & POS38 >= 25.5e6 & POS38 <= 34e6)]
  if (!nrow(s)) return(GRanges())
  reduce(GRanges(as.character(s$chr), IRanges(s$POS38 - 1e6, s$POS38 + 1e6)))
}
nl <- function(d, lp) length(loci(d, lp))
found <- function(d, a, b) { la <- loci(d, a); if (!length(la)) 0L else sum(overlapsAny(la, loci(d, b))) }

pop_lp <- function(x) {
  r <- attr(x, "r_pheno"); n <- attr(x, "n_lab"); N <- attr(x, "n_unlab")
  zz <- function(b, v) b / sqrt(v)
  mm <- cal_qt_opt_vec(data.frame(Z1 = zz(x[, "Z1mm_beta"], x[, "Z1mm_var"]), Z2 = zz(x[, "Z2mm_beta"], x[, "Z2mm_var"]),
                                  Z3 = zz(x[, "ObsG_hat_beta"], x[, "ObsG_var_hat_beta"])),
                       n = n * x[, "Z3lm_var"] / x[, "ObsG_var_hat_beta"],
                       N = N * x[, "Z2lm_var"] / x[, "Z2mm_var"], r = r, eaf = x[, "eaf"])
  lm <- cal_qt_opt_vec(data.frame(Z1 = zz(x[, "Z1lm_beta"], x[, "Z1lm_var"]), Z2 = zz(x[, "Z2lm_beta"], x[, "Z2lm_var"]),
                                  Z3 = zz(x[, "Z3lm_beta"], x[, "Z3lm_var"])), n = n, N = N, r = r, eaf = x[, "eaf"])
  data.table(rsid = rownames(x), SynPALM = x[, "SynSurrG_negative_log10_pvalues"], ObsGWAS = x[, "ObsG_negative_log10_pvalues"],
             POP_mm = z_to_lp(mm$Z), POP_lm = z_to_lp(lm$Z), ObsLM = z_to_lp(zz(x[, "Z3lm_beta"], x[, "Z3lm_var"])),
             se_syn = sqrt(x[, "SynSurrG_var_hat_beta"]))
}
M <- c("SynPALM", "ObsGWAS", "POP_mm", "POP_lm", "ObsLM")

pop_rows <- list(); lab_rows <- list(); ss_rows <- list()
for (i in seq_len(nrow(pr))) {
  p <- pr$protein[i]
  w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
  pos <- w[, .(rsid, chr = as.character(chr), POS38)]
  ff <- sprintf("%s%s_%s.rds", OUT, p, c("full", "full_jointPERM1", "labonly"))
  ## (A) POP-GWAS, real data and joint permutation
  for (k in 1:2) if (file.exists(ff[k])) {
    x <- readRDS(ff[k]); d <- merge(pop_lp(x), pos, by = "rsid")
    chk <- if (k == 1) merge(d[, .(rsid, SynPALM, ObsGWAS)], w[, .(rsid, a = SynSurrG_negative_log10_pvalues, b = ObsG_negative_log10_pvalues)], by = "rsid")[
      , max(abs(SynPALM - a), abs(ObsGWAS - b))] else NA_real_
    row <- data.table(protein = p, q = pr$q[i], data = c("real", "jointPERM1")[k], rho_oof = attr(x, "rho_oof"),
                      r_pheno = attr(x, "r_pheno"), n_lab = attr(x, "n_lab"), n_unlab = attr(x, "n_unlab"),
                      n_snps = nrow(d), check_vs_main = chk,
                      ## relative efficiency: median chi-square ratio at the SNPs genome-wide significant
                      ## for Syn-PALM or POP_mm (z-scale, so free of the per-allele SE scaling of Fun.R)
                      chisq_ratio_POPmm_SynPALM = d[SynPALM >= THR | POP_mm >= THR, median(lp_to_chi(POP_mm) / lp_to_chi(SynPALM))],
                      chisq_ratio_ObsGWAS_SynPALM = d[SynPALM >= THR | ObsGWAS >= THR, median(lp_to_chi(ObsGWAS) / lp_to_chi(SynPALM))])
    for (m in M) {
      row[[paste0("lambda_", m)]] <- d[rsid %in% rnd, lam(get(m))]
      row[[paste0("gws_", m)]]    <- sum(d[[m]] >= THR, na.rm = TRUE)
      row[[paste0("loci_", m)]]   <- nl(d, m)
    }
    row$loci_SynPALM_in_POPmm <- found(d, "SynPALM", "POP_mm")
    row$loci_POPmm_in_SynPALM <- found(d, "POP_mm", "SynPALM")
    pop_rows[[length(pop_rows) + 1]] <- row
  }
  ## (B) labelled-only Syn-PALM against the full (main) analysis
  if (file.exists(ff[3])) {
    x <- readRDS(ff[3])
    d <- merge(w[, .(rsid, chr = as.character(chr), POS38, full_lp = SynSurrG_negative_log10_pvalues,
                     full_se = sqrt(SynSurrG_var_hat_beta), obs_lp_main = ObsG_negative_log10_pvalues,
                     obs_b_main = ObsG_hat_beta, obs_se_main = sqrt(ObsG_var_hat_beta))],
               data.table(rsid = rownames(x), lab_lp = x[, "SynSurrG_negative_log10_pvalues"],
                          lab_se = sqrt(x[, "SynSurrG_var_hat_beta"]), obs_lp = x[, "ObsG_negative_log10_pvalues"],
                          obs_b = x[, "ObsG_hat_beta"], obs_se = sqrt(x[, "ObsG_var_hat_beta"])), by = "rsid")
    lab_rows[[length(lab_rows) + 1]] <- data.table(
      protein = p, q = pr$q[i], rho_oof = attr(x, "rho_oof"), n_lab = attr(x, "n_lab"), n_snps = nrow(d),
      obs_check_maxdiff_lp = d[, max(abs(obs_lp - obs_lp_main))], obs_check_maxreldiff_se = d[, max(abs(obs_se / obs_se_main - 1))],
      se_ratio_full = d[rsid %in% rnd, median(full_se / obs_se_main)], se_ratio_lab = d[rsid %in% rnd, median(lab_se / obs_se)],
      lambda_full = d[rsid %in% rnd, lam(full_lp)], lambda_lab = d[rsid %in% rnd, lam(lab_lp)], lambda_obs = d[rsid %in% rnd, lam(obs_lp)],
      gws_obs = sum(d$obs_lp >= THR), gws_lab = sum(d$lab_lp >= THR), gws_full = sum(d$full_lp >= THR),
      loci_obs = nl(d, "obs_lp"), loci_lab = nl(d, "lab_lp"), loci_full = nl(d, "full_lp"),
      tau_TS_lab = attr(x, "vc")$used[["tau_TS"]], sigma_TS_lab = attr(x, "vc")$used[["sigma_TS"]])
  }
  ## SynSurr (unrelated subset) from the main genome-wide scan
  sc <- readRDS(sprintf("results/allukb_protein_pQTL_european_discovery/%s/%s_cf_noPRS_noPC_scan_allSNPs.rds", p, p))
  sc <- sc[rownames(sc) %in% pos$rsid, , drop = FALSE]   # autosomal SNPs of the main results
  d <- data.table(rsid = rownames(sc), SynPALM = sc[, "SynSurrG_negative_log10_pvalues"], ObsGWAS = sc[, "ObsG_negative_log10_pvalues"],
                  SynSurr_unrel = sc[, "SynSurr_negative_log10_pvalues"], Obs_unrel = sc[, "Obs_negative_log10_pvalues"],
                  se_ratio = sqrt(sc[, "SynSurr_var_hat_beta"] / sc[, "SynSurrG_var_hat_beta"]))
  d <- merge(d, pos, by = "rsid")
  row <- data.table(protein = p, q = pr$q[i], n_snps = nrow(d), se_ratio_SynSurr_SynPALM = d[rsid %in% rnd, median(se_ratio, na.rm = TRUE)])
  for (m in c("SynPALM", "ObsGWAS", "SynSurr_unrel", "Obs_unrel")) {
    row[[paste0("gws_", m)]] <- sum(d[[m]] >= THR, na.rm = TRUE); row[[paste0("loci_", m)]] <- nl(d, m)
  }
  ss_rows[[length(ss_rows) + 1]] <- row
  cat(p, "\n")
}
pop <- rbindlist(pop_rows, fill = TRUE); lab <- rbindlist(lab_rows); ss <- rbindlist(ss_rows)
fwrite(pop, paste0(OUT, "pop_per_protein.csv")); fwrite(lab, paste0(OUT, "labonly_per_protein.csv"))
fwrite(ss, paste0(OUT, "synsurr_unrelated_per_protein.csv"))
if (nrow(pop)) {
  ps <- pop[, c(list(proteins = .N, check_vs_main = max(check_vs_main), r_pheno_median = median(r_pheno),
                     chisq_ratio_POPmm_SynPALM = median(chisq_ratio_POPmm_SynPALM, na.rm = TRUE),
                     chisq_ratio_ObsGWAS_SynPALM = median(chisq_ratio_ObsGWAS_SynPALM, na.rm = TRUE)),
                setNames(lapply(M, function(m) median(get(paste0("lambda_", m)))), paste0("lambda_med_", M)),
                setNames(lapply(M, function(m) max(get(paste0("lambda_", m)))), paste0("lambda_max_", M)),
                setNames(lapply(M, function(m) sum(get(paste0("gws_", m)))), paste0("gws_", M)),
                setNames(lapply(M, function(m) sum(get(paste0("loci_", m)))), paste0("loci_", M)),
                list(loci_SynPALM_in_POPmm = sum(loci_SynPALM_in_POPmm), loci_POPmm_in_SynPALM = sum(loci_POPmm_in_SynPALM))), by = data]
  fwrite(ps, paste0(OUT, "pop_summary.csv")); print(t(ps))
}
if (nrow(lab)) {
  ls <- lab[, .(proteins = .N, obs_check_lp = max(obs_check_maxdiff_lp), obs_check_se = max(obs_check_maxreldiff_se),
                se_ratio_full_med = median(se_ratio_full), se_ratio_lab_med = median(se_ratio_lab),
                gain_share_from_labelled = median((1 - se_ratio_lab) / (1 - se_ratio_full)),
                gws_obs = sum(gws_obs), gws_lab = sum(gws_lab), gws_full = sum(gws_full),
                loci_obs = sum(loci_obs), loci_lab = sum(loci_lab), loci_full = sum(loci_full))]
  fwrite(ls, paste0(OUT, "labonly_summary.csv")); print(t(ls))
}
print(ss[, lapply(.SD, sum), .SDcols = patterns("^gws|^loci")]); print(median(ss$se_ratio_SynSurr_SynPALM))

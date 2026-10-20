## Summary of the REML refit (code/post_cf3/vcrefit/vcrefit.R) against the main (MoM) run.
## For each protein of results/vcrefit_cf_noPRS_noPC/proteins.csv, on the protein's SNP set
## (main P < 1e-3 for Syn-PALM or Obs-GWAS + 20,000 random SNPs; autosomes):
##   MoM vs REML variance components (raw; used after the PSD constraint; projection flags);
##   correlation of -log10 P (Syn-PALM; Obs-GWAS) over the 20,000 random SNPs and over all
##   scanned SNPs; lambda_GC over the random SNPs; genome-wide significant SNPs (P < 5e-8,
##   MHC chr6:25.5-34 Mb removed) and primary associations (+-1 Mb merged windows) of the
##   main and REML runs and their overlap; the change of the estimate at the main GWS SNPs
##   in units of its main SE and the ratio of the SEs.
## main = results/web_SynGWAS_summary_results_cf_noPRS_noPC/<p>_negative_log10_pvals_v2geno_crossfit.rds
##   Rscript code/post_cf3/vcrefit/summarize_vcrefit.R
## Writes results/vcrefit_cf_noPRS_noPC/{vc_mom_vs_reml.csv, vcrefit_per_protein.csv, vcrefit_summary.csv}
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
setwd("<PROJECT_DIR>")
OUT <- "results/vcrefit_cf_noPRS_noPC/"
pr  <- fread(paste0(OUT, "proteins.csv"))
rnd <- readLines("results/sens_cf_noPRS_noPC/snps_random20k.txt")
THR <- -log10(5e-8); q50 <- qchisq(0.5, 1)
lam <- function(lp) qchisq(10^(-median(lp)), 1, lower.tail = FALSE) / q50
loci <- function(d, lp) {
  s <- d[get(lp) >= THR & !(chr == 6 & POS38 >= 25.5e6 & POS38 <= 34e6)]
  if (!nrow(s)) return(GRanges())
  reduce(GRanges(as.character(s$chr), IRanges(s$POS38 - 1e6, s$POS38 + 1e6)))
}
nm <- c("tau_T2", "tau_TS", "tau_S2", "sigma_T2", "sigma_TS", "sigma_S2")
vcrows <- list(); rows <- list(); miss <- character()
for (i in seq_len(nrow(pr))) {
  p <- pr$protein[i]
  vf <- sprintf("%svc/%s_reml.rds", OUT, p); sf <- sprintf("%sscan/%s_reml_scan_subset.rds", OUT, p)
  if (!file.exists(vf) || !file.exists(sf)) { miss <- c(miss, p); next }
  V <- readRDS(vf); x <- readRDS(sf); vcr <- attr(x, "vc")
  L <- readRDS(sprintf("results/table_results_cf_noPRS_noPC/%s_loci.rds", p)); vcm <- L$vc
  stopifnot(max(abs(V$mom_raw[nm] - vcm$raw[nm])) < 1e-8)       # MoM recomputed = production
  vcrows[[length(vcrows) + 1]] <- data.table(protein = p, set = pr$set[i], rho_oof = pr$rho_oof[i],
    vc = nm, mom_raw = vcm$raw[nm], mom_used = vcm$used[nm], reml_raw = V$reml_raw[nm], reml_used = vcr$used[nm])
  w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
  d <- merge(w[, .(rsid, chr, POS38, m_lp = SynSurrG_negative_log10_pvalues, m_b = SynSurrG_hat_beta,
                   m_se = sqrt(SynSurrG_var_hat_beta), mo_lp = ObsG_negative_log10_pvalues,
                   mo_b = ObsG_hat_beta, mo_se = sqrt(ObsG_var_hat_beta))],
             data.table(rsid = rownames(x), n_lp = x[, "SynSurrG_negative_log10_pvalues"],
                        n_b = x[, "SynSurrG_hat_beta"], n_se = sqrt(x[, "SynSurrG_var_hat_beta"]),
                        no_lp = x[, "ObsG_negative_log10_pvalues"], no_b = x[, "ObsG_hat_beta"],
                        no_se = sqrt(x[, "ObsG_var_hat_beta"])), by = "rsid")
  stopifnot(nrow(d) == nrow(x), all(rnd %in% d$rsid))
  ## underflow guard as post_protein.R
  for (k in c("n", "no")) { lp <- paste0(k, "_lp"); bad <- which(!is.finite(d[[lp]]))
    if (length(bad)) set(d, bad, lp, -pchisq((d[[paste0(k, "_b")]][bad] / d[[paste0(k, "_se")]][bad])^2, 1,
                                              lower.tail = FALSE, log.p = TRUE) / log(10)) }
  dr <- d[rsid %in% rnd]
  nomhc <- d[!(chr == 6 & POS38 >= 25.5e6 & POS38 <= 34e6)]
  gm <- nomhc[m_lp >= THR, rsid]; gn <- nomhc[n_lp >= THR, rsid]
  go <- nomhc[mo_lp >= THR, rsid]; gon <- nomhc[no_lp >= THR, rsid]
  lm_ <- loci(d, "m_lp"); ln <- loci(d, "n_lp"); lo <- loci(d, "mo_lp"); lon <- loci(d, "no_lp")
  gws <- d[m_lp >= THR]
  rows[[length(rows) + 1]] <- data.table(
    protein = p, set = pr$set[i], q = pr$q[i], rho_oof = pr$rho_oof[i],
    raw_tauT2_le0 = vcm$raw[["tau_T2"]] <= 0,
    mom_projected_G = vcm$projected_G, reml_projected_G = vcr$projected_G, reml_projected_E = vcr$projected_E,
    reml_converged = all(V$converged),
    mom_rho_g = vcm$rho_g, reml_rho_g = vcr$rho_g,
    mom_tauT2_raw = vcm$raw[["tau_T2"]], reml_tauT2 = V$reml_raw[["tau_T2"]], reml_tauT2_uni = V$reml_uni_T[["tau"]],
    mom_tauTS_raw = vcm$raw[["tau_TS"]], reml_tauTS = V$reml_raw[["tau_TS"]],
    mom_tauS2_raw = vcm$raw[["tau_S2"]], reml_tauS2 = V$reml_raw[["tau_S2"]],
    reml_biv_rho_g = V$reml_biv$G[1, 2] / sqrt(max(V$reml_biv$G[1, 1] * V$reml_biv$G[2, 2], 1e-300)),
    ## unconstrained bivariate REML outside the PSD cone? and its -2 log-lik gain over the PSD fit
    free_nonPSD_G = min(eigen(V$reml_biv_free$G, symmetric = TRUE, only.values = TRUE)$values) < 0,
    free_tauT2 = V$reml_biv_free$G[1, 1],
    free_rho_g = V$reml_biv_free$G[1, 2] / sqrt(abs(V$reml_biv_free$G[1, 1] * V$reml_biv_free$G[2, 2])),
    free_2ll_gain = 2 * (V$reml_biv$nll - V$reml_biv_free$nll),
    reml_shift_G = max(abs(vcr$used[c("tau_T2", "tau_TS", "tau_S2")] - V$reml_raw[c("tau_T2", "tau_TS", "tau_S2")])),
    n_snps = nrow(d),
    r_syn_random = cor(dr$m_lp, dr$n_lp), r_syn_all = cor(d$m_lp, d$n_lp),
    r_obs_random = cor(dr$mo_lp, dr$no_lp),
    lambda_syn_main = lam(dr$m_lp), lambda_syn_reml = lam(dr$n_lp),
    lambda_obs_main = lam(dr$mo_lp), lambda_obs_reml = lam(dr$no_lp),
    gws_main = length(gm), gws_reml = length(gn), gws_both = length(intersect(gm, gn)),
    obs_gws_main = length(go), obs_gws_reml = length(gon),
    loci_main = length(lm_), loci_reml = length(ln),
    loci_main_found = if (length(lm_)) sum(overlapsAny(lm_, ln)) else 0L,
    loci_reml_found = if (length(ln)) sum(overlapsAny(ln, lm_)) else 0L,
    obs_loci_main = length(lo), obs_loci_reml = length(lon),
    dz_med = if (nrow(gws)) gws[, median(abs(n_b - m_b) / m_se)] else NA_real_,
    dz_max = if (nrow(gws)) gws[, max(abs(n_b - m_b) / m_se)] else NA_real_,
    se_ratio_med = if (nrow(gws)) gws[, median(n_se / m_se)] else NA_real_,
    beta_ratio_med = if (nrow(gws)) gws[, median(n_b / m_b)] else NA_real_,
    max_abs_dlp = d[, max(abs(n_lp - m_lp))])
}
if (length(miss)) cat("MISSING outputs for", length(miss), "proteins:", paste(miss, collapse = " "), "\n")
vcd <- rbindlist(vcrows); r <- rbindlist(rows)
fwrite(vcd, paste0(OUT, "vc_mom_vs_reml.csv")); fwrite(r, paste0(OUT, "vcrefit_per_protein.csv"))
agg <- function(x) x[, .(proteins = .N, reml_projected_G = sum(reml_projected_G), reml_shift_G_median = median(reml_shift_G),
  reml_shift_G_max = max(reml_shift_G), reml_not_converged = sum(!reml_converged),
  free_reml_nonPSD = sum(free_nonPSD_G), free_2ll_gain_median = median(free_2ll_gain), free_2ll_gain_max = max(free_2ll_gain),
  r_syn_random_median = median(r_syn_random), r_syn_random_min = min(r_syn_random),
  r_syn_all_median = median(r_syn_all), r_syn_all_min = min(r_syn_all),
  lambda_syn_main_median = median(lambda_syn_main), lambda_syn_reml_median = median(lambda_syn_reml),
  lambda_syn_absdiff_max = max(abs(lambda_syn_reml - lambda_syn_main)),
  lambda_obs_main_median = median(lambda_obs_main), lambda_obs_reml_median = median(lambda_obs_reml),
  gws_main = sum(gws_main), gws_reml = sum(gws_reml), gws_overlap_of_main = sum(gws_both) / sum(gws_main),
  gws_overlap_of_reml = sum(gws_both) / sum(gws_reml),
  loci_main = sum(loci_main), loci_reml = sum(loci_reml),
  loci_main_recovered = sum(loci_main_found) / sum(loci_main), loci_reml_in_main = sum(loci_reml_found) / sum(loci_reml),
  obs_loci_main = sum(obs_loci_main), obs_loci_reml = sum(obs_loci_reml),
  syn_obs_ratio_main = sum(loci_main) / sum(obs_loci_main), syn_obs_ratio_reml = sum(loci_reml) / sum(obs_loci_reml),
  dz_median = median(dz_med, na.rm = TRUE), dz_max = max(dz_max, na.rm = TRUE),
  se_ratio_median = median(se_ratio_med, na.rm = TRUE), beta_ratio_median = median(beta_ratio_med, na.rm = TRUE))]
sm <- rbind(cbind(group = "projected", agg(r[set == "projected"])),
            cbind(group = "projected, raw tau_T2 <= 0", agg(r[set == "projected" & raw_tauT2_le0 == TRUE])),
            cbind(group = "projected, raw tau_T2 > 0", agg(r[set == "projected" & raw_tauT2_le0 == FALSE])),
            cbind(group = "reference (not projected)", agg(r[set == "reference"])))
fwrite(sm, paste0(OUT, "vcrefit_summary.csv"))
vs <- vcd[, .(median_mom_raw = median(mom_raw), median_mom_used = median(mom_used), median_reml = median(reml_raw),
              median_abs_diff_used = median(abs(reml_used - mom_used)), r_used = cor(reml_used, mom_used)), by = .(set, vc)]
fwrite(vs, paste0(OUT, "vc_mom_vs_reml_summary.csv"))
options(width = 250)
print(r, digits = 3); print(t(sm)); print(vs, digits = 3)

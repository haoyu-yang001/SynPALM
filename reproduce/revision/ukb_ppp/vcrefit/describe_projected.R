## Revision (NC R1): are the method-of-moments (MoM) variance components adequate?
## Descriptive part. For every protein of the main analysis (cross-fitted, no PRS,
## no PCs) collect, from results/table_results_cf_noPRS_noPC/<p>_loci.rds,
## the raw and used MoM variance components and the projection flags, rho_oof and
## the primary associations (P < 5e-8, +-1 Mb merged, MHC removed) of Syn-PALM and
## Obs-GWAS; and, from the web summary statistics
## results/web_SynGWAS_summary_results_cf_noPRS_noPC/<p>_negative_log10_pvals_v2geno_crossfit.rds
## (autosomes only), the genomic inflation factor lambda_GC = median chi2 / qchisq(0.5, 1)
## of both methods over all autosomal SNPs and over the 20,000 random SNPs of
## results/sens_cf_noPRS_noPC/snps_random20k.txt. Then tabulate projected vs not projected.
##   Rscript code/post_cf3/vcrefit/describe_projected.R       (SLURM: submit_describe.sh)
## Writes results/vcrefit_cf_noPRS_noPC/describe_per_protein.csv and describe_summary.csv
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
OUT <- "results/vcrefit_cf_noPRS_noPC/"; dir.create(OUT, showWarnings = FALSE)
NC  <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
per <- fread("results/table_results_cf_noPRS_noPC/summary_per_protein.csv")
rnd <- readLines("results/sens_cf_noPRS_noPC/snps_random20k.txt")
q50 <- qchisq(0.5, 1)
lam <- function(lp) qchisq(10^(-median(lp, na.rm = TRUE)), 1, lower.tail = FALSE) / q50
pp_file <- paste0(OUT, "describe_per_protein.csv")
if (!file.exists(pp_file)) {
  ## proteins are processed in chunks of 200 saved in tmp_describe/, so a preempted job resumes
  TMPD <- paste0(OUT, "tmp_describe/"); dir.create(TMPD, showWarnings = FALSE)
  one <- function(p) {
    L <- readRDS(sprintf("results/table_results_cf_noPRS_noPC/%s_loci.rds", p))
    w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
    stopifnot(all(w$chr %in% 1:22))
    wr <- w[rsid %in% rnd]
    v  <- L$vc
    as.data.table(c(list(protein = p, rho_oof = L$rho_oof, n_obs = L$n_obs,
               projected_G = v$projected_G, projected_E = v$projected_E,
               truncated_G = v$detail_G$truncated, shrunk_G = v$detail_G$shrunk),
               as.list(setNames(v$raw, paste0("raw_", names(v$raw)))),
               as.list(setNames(v$used, paste0("used_", names(v$used)))),
               list(rho_g_used = v$rho_g, rho_e_used = v$rho_e,
               syn_loci = nrow(L$by_threshold[["5e-08"]]$SynG$loci),
               obs_loci = nrow(L$by_threshold[["5e-08"]]$ObsG$loci),
               syn_snps = L$by_threshold[["5e-08"]]$SynG$n_snps,
               obs_snps = L$by_threshold[["5e-08"]]$ObsG$n_snps,
               n_snps_web = nrow(w), n_random = nrow(wr),
               lambda_syn = lam(w$SynSurrG_negative_log10_pvalues),
               lambda_obs = lam(w$ObsG_negative_log10_pvalues),
               lambda_syn_random = lam(wr$SynSurrG_negative_log10_pvalues),
               lambda_obs_random = lam(wr$ObsG_negative_log10_pvalues))))
  }
  chunks <- split(per$protein, ceiling(seq_along(per$protein) / 200))
  for (k in seq_along(chunks)) {
    cf <- sprintf("%schunk%02d.rds", TMPD, k)
    if (file.exists(cf)) next
    rows <- parallel::mclapply(chunks[[k]], function(p) tryCatch(one(p), error = function(e) e),
                               mc.cores = NC, mc.preschedule = FALSE)
    bad <- vapply(rows, function(x) !is.data.table(x), logical(1))
    if (any(bad)) stop(sum(bad), " proteins failed, first: ", chunks[[k]][which(bad)[1]], ": ",
                       if (inherits(rows[[which(bad)[1]]], "error")) conditionMessage(rows[[which(bad)[1]]]) else "worker killed?")
    saveRDS(rbindlist(rows), cf); cat("chunk", k, "/", length(chunks), "done\n")
  }
  d <- rbindlist(lapply(seq_along(chunks), function(k) readRDS(sprintf("%schunk%02d.rds", TMPD, k))))
  ## raw 2x2 genetic block: correlation (NA when a diagonal is <= 0)
  d[, rho_g_raw := fifelse(raw_tau_T2 > 0 & raw_tau_S2 > 0, raw_tau_TS / sqrt(pmax(raw_tau_T2 * raw_tau_S2, 0)), NA_real_)]
  fwrite(d, pp_file)
}
d <- fread(pp_file)
stopifnot(nrow(d) == nrow(per), !anyNA(d$projected_G))
## sanity checks against summary_per_protein.csv and vc_check
chk <- merge(d, per[, .(protein, s_syn = `SynG_loci_5e-08`, s_obs = `ObsG_loci_5e-08`, s_rho = rho_oof)], by = "protein")
stopifnot(all(chk$syn_loci == chk$s_syn), all(chk$obs_loci == chk$s_obs), max(abs(chk$rho_oof - chk$s_rho)) < 1e-12)
vcc <- fread("results/vc_check_cf_noPRS_noPC.csv")[set == "main" & done == TRUE]
m <- merge(d, vcc[, .(protein, vproj = projected_G, vt = raw_tau_T2)], by = "protein")
cat("vc_check agreement: projected_G", sum(m$projected_G == m$vproj), "/", nrow(m),
    "; max |raw tau_T2 diff|", signif(max(abs(m$raw_tau_T2 - m$vt)), 3), "\n")

d[, group := fifelse(!projected_G, "not projected",
              fifelse(raw_tau_T2 <= 0, "projected, raw tau_T2 <= 0", "projected, raw tau_T2 > 0"))]
summ <- function(x) {
  x[, .(proteins = .N,
        rho_median = median(rho_oof), rho_q1 = quantile(rho_oof, 0.25), rho_q3 = quantile(rho_oof, 0.75),
        rho_lt_0.1 = mean(rho_oof < 0.1),
        tauT2_raw_median = median(raw_tau_T2), tauS2_raw_median = median(raw_tau_S2),
        syn_loci = sum(syn_loci), obs_loci = sum(obs_loci), ratio_loci = sum(syn_loci) / sum(obs_loci),
        prot_with_assoc = sum(syn_loci + obs_loci > 0),
        ratio_loci_median_proteins_with_assoc = median((syn_loci / obs_loci)[obs_loci > 0]),
        frac_more_syn = mean(syn_loci > obs_loci), frac_fewer_syn = mean(syn_loci < obs_loci),
        lambda_syn_median = median(lambda_syn), lambda_syn_q1 = quantile(lambda_syn, 0.25), lambda_syn_q3 = quantile(lambda_syn, 0.75),
        lambda_obs_median = median(lambda_obs), lambda_obs_q1 = quantile(lambda_obs, 0.25), lambda_obs_q3 = quantile(lambda_obs, 0.75),
        lambda_syn_gt_1.1 = mean(lambda_syn > 1.1), lambda_obs_gt_1.1 = mean(lambda_obs > 1.1),
        lambda_syn_lt_0.95 = mean(lambda_syn < 0.95), lambda_obs_lt_0.95 = mean(lambda_obs < 0.95),
        lambda_ratio_median = median(lambda_syn / lambda_obs))]
}
tab <- rbind(cbind(group = "all", summ(d)),
             cbind(group = "not projected", summ(d[projected_G == FALSE])),
             cbind(group = "projected (all)", summ(d[projected_G == TRUE])),
             cbind(group = "projected, raw tau_T2 <= 0", summ(d[projected_G & raw_tau_T2 <= 0])),
             cbind(group = "projected, raw tau_T2 > 0", summ(d[projected_G & raw_tau_T2 > 0])))
## the same within rho_oof strata (projection is more frequent at low rho?)
d[, rho_bin := cut(rho_oof, c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, Inf), right = FALSE,
                   labels = c("<0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", ">=0.5"))]
strat <- d[, .(proteins = .N, frac_projected = mean(projected_G),
               ratio_loci_proj = sum(syn_loci[projected_G]) / sum(obs_loci[projected_G]),
               ratio_loci_notproj = sum(syn_loci[!projected_G]) / sum(obs_loci[!projected_G]),
               lambda_syn_proj = median(lambda_syn[projected_G]), lambda_syn_notproj = median(lambda_syn[!projected_G]),
               lambda_obs_proj = median(lambda_obs[projected_G]), lambda_obs_notproj = median(lambda_obs[!projected_G])),
           by = rho_bin][order(rho_bin)]
## projection mode: which constraint was hit
modes <- d[projected_G == TRUE, .N, by = .(raw_tauT2_le0 = raw_tau_T2 <= 0, abs_rho_g_raw_gt1 = !is.na(rho_g_raw) & abs(rho_g_raw) > 1)]
wt <- d[, .(w_lambda_syn = wilcox.test(lambda_syn[projected_G], lambda_syn[!projected_G])$p.value,
            w_lambda_obs = wilcox.test(lambda_obs[projected_G], lambda_obs[!projected_G])$p.value,
            w_rho = wilcox.test(rho_oof[projected_G], rho_oof[!projected_G])$p.value)]
fwrite(tab, paste0(OUT, "describe_summary.csv"))
fwrite(strat, paste0(OUT, "describe_by_rho.csv"))
options(width = 250)
print(tab, digits = 3); print(strat, digits = 3); print(modes); print(wt, digits = 3)

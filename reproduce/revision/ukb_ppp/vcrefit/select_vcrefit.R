## Revision (NC R1): REML refit of the Syn-PALM variance components -- protein and SNP selection.
## 50 proteins whose MoM genetic block was projected onto the PSD cone in the main
## analysis (projected_G), drawn at random (seed 20261020) stratified by quintile of the
## out-of-fold correlation rho_oof among the projected proteins and by the sign of the raw
## MoM tau_T2 (5 per quintile x {raw tau_T2 <= 0, > 0}; a stratum with fewer than 5 is
## filled from the other sign of the same quintile), plus 10 non-projected reference
## proteins (2 per quintile of rho_oof among the non-projected). Projection flags and raw
## VCs from results/vcrefit_cf_noPRS_noPC/describe_per_protein.csv (= the vc stored with
## each protein's main scan).
## SNP set per protein (as code/post_cf3/sens/select_sens.R): autosomal SNPs with P < 1e-3
## for Syn-PALM or Obs-GWAS in the main analysis, plus the 20,000 random autosomal SNPs of
## results/sens_cf_noPRS_noPC/snps_random20k.txt.
##   Rscript code/post_cf3/vcrefit/select_vcrefit.R
## Writes results/vcrefit_cf_noPRS_noPC/{proteins.csv, snps/snps_<protein>.txt}
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
OUT <- "results/vcrefit_cf_noPRS_noPC/"; dir.create(paste0(OUT, "snps"), recursive = TRUE, showWarnings = FALSE)
d <- fread(paste0(OUT, "describe_per_protein.csv"))
set.seed(20261020)
pj <- d[projected_G == TRUE]
pj[, q := cut(rho_oof, quantile(rho_oof, 0:5 / 5), include.lowest = TRUE, labels = 1:5)]
pj[, neg := raw_tau_T2 <= 0]
sel <- rbindlist(lapply(1:5, function(k) {
  a <- pj[q == k & neg]; b <- pj[q == k & !neg]
  na <- min(5L, nrow(a)); nb <- min(5L, nrow(b))
  if (na < 5) nb <- min(10L - na, nrow(b)); if (nb < 5) na <- min(10L - nb, nrow(a))
  rbind(a[sample(.N, na)], b[sample(.N, nb)])
}))
np <- d[projected_G == FALSE]
np[, q := cut(rho_oof, quantile(rho_oof, 0:5 / 5), include.lowest = TRUE, labels = 1:5)]
np[, neg := raw_tau_T2 <= 0]
ref <- np[, .SD[sample(.N, 2)], by = q]
sel <- rbind(cbind(set = "projected", sel), cbind(set = "reference", ref), fill = TRUE)
hdr <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
sel[, index := match(protein, hdr)]
if (anyNA(sel$index)) sel[is.na(index), index := match(protein, gsub("-", "_", tolower(hdr)))]
stopifnot(!anyNA(sel$index), !anyDuplicated(sel$protein))
sel <- sel[order(set, q, -rho_oof)]
fwrite(sel[, .(set, q, protein, index, rho_oof, raw_tau_T2, raw_tau_TS, raw_tau_S2, rho_g_raw,
               syn_loci, obs_loci, lambda_syn, lambda_obs)], paste0(OUT, "proteins.csv"))
print(sel[, .(set, q, protein, index, rho_oof, raw_tau_T2, syn_loci, obs_loci)])
rnd <- readLines("results/sens_cf_noPRS_noPC/snps_random20k.txt")
for (p in sel$protein) {
  w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
  stopifnot(all(rnd %in% w$rsid))
  cand <- w[SynSurrG_negative_log10_pvalues >= 3 | ObsG_negative_log10_pvalues >= 3, rsid]
  writeLines(union(cand, rnd), sprintf("%ssnps/snps_%s.txt", OUT, p))
  cat(p, length(cand), "candidate SNPs\n")
}

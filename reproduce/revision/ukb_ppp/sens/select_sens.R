## Proteins and SNP sets of the sensitivity analysis (revision 2026-10-19; reply to
## Reviewer 1, Comments 4 and 8). Six proteins are drawn at random (seed 20261019)
## from each quartile of the out-of-fold correlation of the main analysis, among the
## proteins with at least one Syn-PALM or Obs-GWAS primary association. For each, the
## scanned SNPs are all autosomal SNPs with P < 1e-3 for Syn-PALM or Obs-GWAS in the
## main analysis, plus 20,000 random autosomal SNPs (the same for every protein) for
## the genome-wide correlation of -log10 P.
##   Rscript code/post_cf3/sens/select_sens.R
## Writes results/sens_cf_noPRS_noPC/{proteins.csv, snps_<protein>.txt}
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
OUT <- "results/sens_cf_noPRS_noPC/"; dir.create(OUT, showWarnings = FALSE)
per <- fread("results/table_results_cf_noPRS_noPC/summary_per_protein.csv")
per <- per[`SynG_loci_5e-08` + `ObsG_loci_5e-08` > 0]
per[, q := cut(rho_oof, quantile(rho_oof, 0:4 / 4), include.lowest = TRUE, labels = 1:4)]
set.seed(20261019)
sel <- per[, .SD[sample(.N, 6)], by = q][order(q, -rho_oof)]
## protein index (column position) used by the step scripts
hdr <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
sel[, index := match(protein, gsub("-", "_", tolower(hdr)))]
if (anyNA(sel$index)) sel[is.na(index), index := match(protein, tolower(hdr))]
stopifnot(!anyNA(sel$index))
fwrite(sel[, .(q, protein, index, rho_oof, SynG = `SynG_loci_5e-08`, ObsG = `ObsG_loci_5e-08`)], paste0(OUT, "proteins.csv"))
print(sel[, .(q, protein, index, rho_oof, `SynG_loci_5e-08`, `ObsG_loci_5e-08`)])
web1 <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", sel$protein[1]))
set.seed(1); rnd <- sample(web1$rsid, 20000)
for (p in sel$protein) {
  w <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", p))
  cand <- w[SynSurrG_negative_log10_pvalues >= 3 | ObsG_negative_log10_pvalues >= 3, rsid]
  writeLines(union(cand, rnd), sprintf("%ssnps_%s.txt", OUT, p))
  cat(p, length(cand), "candidate SNPs\n")
}
writeLines(rnd, paste0(OUT, "snps_random20k.txt"))

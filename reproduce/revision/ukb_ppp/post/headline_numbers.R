## Numbers quoted in the Results of the revised manuscript (UKB-PPP, no PRS /
## no PCs), all at P < 5e-8, autosomes, from the files of aggregate_all.R.
##   Rscript code/post_cf3/headline_numbers.R
## Writes results/table_results_cf_noPRS_noPC/headline_numbers.txt
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
D <- "results/table_results_cf_noPRS_noPC/"
sink(paste0(D, "headline_numbers.txt"), split = TRUE)
per <- fread(paste0(D, "summary_per_protein.csv")); setnames(per, gsub("-", "", names(per)))
f <- function(x) format(x, big.mark = ",")
cat("proteins:", nrow(per), "\n")
cat(sprintf("out-of-fold rho: median %.3f, quartiles %.3f %.3f, 10%% %.3f, 90%% %.3f\n", median(per$rho_oof),
            quantile(per$rho_oof, .25), quantile(per$rho_oof, .75), quantile(per$rho_oof, .1), quantile(per$rho_oof, .9)))
cat(sprintf("n_obs: median %s (range %s-%s)\n", f(median(per$n_obs)), f(min(per$n_obs)), f(max(per$n_obs))))
s <- function(col) sum(per[[col]], na.rm = TRUE)
for (k in c("snps_noMHC", "snps_all", "loci"))
  cat(sprintf("%-11s Syn-PALM %s  Obs-GWAS %s  ratio %.2f\n", k, f(s(paste0("SynG_", k, "_5e08"))),
              f(s(paste0("ObsG_", k, "_5e08"))), s(paste0("SynG_", k, "_5e08")) / s(paste0("ObsG_", k, "_5e08"))))
for (m in c("SynG", "ObsG")) {
  z <- fread(paste0(D, "loci_", m, "_5e-08.csv"))
  cat(sprintf("%s primary associations %s: cis %s, trans %s; cis within 50 kb / 100 kb of the gene: %.1f%% / %.1f%%; proteins with >= 1: %s\n",
              m, f(nrow(z)), f(sum(z$cis_trans == "cis", na.rm = TRUE)), f(sum(z$cis_trans == "trans", na.rm = TRUE)),
              100 * z[cis_trans == "cis", mean(abs(POS38 - gene_pos) <= 5e4)],
              100 * z[cis_trans == "cis", mean(abs(POS38 - gene_pos) <= 1e5)], f(uniqueN(z$protein))))
}
top <- per[order(-rho_oof)][1:20]
cat(sprintf("top-20 proteins by rho (%.3f-%.3f): mean primary associations Syn-PALM %.2f, Obs-GWAS %.2f (ratio %.2f); mean significant SNPs (no MHC) %.2f vs %.2f (ratio %.2f)\n",
            min(top$rho_oof), max(top$rho_oof), mean(top$SynG_loci_5e08), mean(top$ObsG_loci_5e08),
            sum(top$SynG_loci_5e08) / sum(top$ObsG_loci_5e08), mean(top$SynG_snps_noMHC_5e08),
            mean(top$ObsG_snps_noMHC_5e08), sum(top$SynG_snps_noMHC_5e08) / sum(top$ObsG_snps_noMHC_5e08)))
print(top[, .(protein, rho = round(rho_oof, 3), n_obs, SynG_loci = SynG_loci_5e08, ObsG_loci = ObsG_loci_5e08,
              SynG_snps = SynG_snps_noMHC_5e08, ObsG_snps = ObsG_snps_noMHC_5e08)])
cat(sprintf("proteins with more / fewer / equal primary associations in Syn-PALM: %d / %d / %d\n",
            per[SynG_loci_5e08 > ObsG_loci_5e08, .N], per[SynG_loci_5e08 < ObsG_loci_5e08, .N], per[SynG_loci_5e08 == ObsG_loci_5e08, .N]))
if (file.exists(paste0(D, "gain_by_rho_strata.csv"))) print(fread(paste0(D, "gain_by_rho_strata.csv")), digits = 3)
sink()

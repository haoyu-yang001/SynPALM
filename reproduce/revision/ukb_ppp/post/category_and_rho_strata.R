## Reviewer 3 comment 2 (gain by prediction accuracy) and minor comment 3
## (Fig 4D/E: category sizes, per-protein rates, enrichment).
##   Rscript code/post_cf3/category_and_rho_strata.R
## Input: results/table_results_cf_noPRS_noPC/summary_per_protein.csv (aggregate_all.R).
## Outputs in results/table_results_cf_noPRS_noPC/:
##   gain_by_rho_strata.csv  proteins, loci and SNPs of Syn-PALM vs Obs-GWAS by rho bin
##   category_rates.csv      per Olink panel: proteins, loci per protein (both methods),
##                           additional loci per protein, and a test of whether the
##                           additional loci are over-represented in the panel
##                           (Poisson GLM of additional loci on panel with the Obs-GWAS
##                           loci as exposure, and a Fisher test of panel vs the rest)
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
D <- "results/table_results_cf_noPRS_noPC/"
per <- fread(paste0(D, "summary_per_protein.csv"))
T <- "5e-08"
per[, `:=`(syn_loci = get(paste0("SynG_loci_", T)), obs_loci = get(paste0("ObsG_loci_", T)),
           syn_snps = get(paste0("SynG_snps_noMHC_", T)), obs_snps = get(paste0("ObsG_snps_noMHC_", T)))]

## ---- by prediction accuracy -------------------------------------------------
br <- c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, Inf)
per[, rho_bin := cut(rho_oof, br, labels = c("<0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", ">=0.6"))]
g <- per[, .(proteins = .N, median_rho = median(rho_oof), syn_loci = sum(syn_loci), obs_loci = sum(obs_loci),
             syn_snps = sum(syn_snps), obs_snps = sum(obs_snps),
             proteins_more_loci = sum(syn_loci > obs_loci), proteins_fewer_loci = sum(syn_loci < obs_loci)), by = rho_bin]
g[, `:=`(loci_ratio = syn_loci / obs_loci, snp_ratio = syn_snps / obs_snps)]
setorder(g, rho_bin); fwrite(g, paste0(D, "gain_by_rho_strata.csv")); print(g, digits = 3)
ct <- cor.test(per$rho_oof, (per$syn_loci + 0.5) / (per$obs_loci + 0.5), method = "spearman", exact = FALSE)
cat(sprintf("Spearman correlation of rho_oof with the per-protein loci ratio: %.3f (P = %.2g)\n", ct$estimate, ct$p.value))

## ---- by functional category (Olink panel, as in Fig 4D) -------------------------
pm <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
pm[, `:=`(key = gsub("-", "_", id), panel4 = sub("_II$", "", panel))]
per[, key := gsub("-", "_", protein)]   # e.g. hla-e in the results, hla_e in the protein map
per[pm, on = .(key = key), panel4 := i.panel4]
per[, extra := pmax(syn_loci - obs_loci, 0)]
cr <- per[!is.na(panel4), .(proteins = .N, syn_loci = sum(syn_loci), obs_loci = sum(obs_loci), extra_loci = sum(extra),
                            syn_loci_per_protein = mean(syn_loci), obs_loci_per_protein = mean(obs_loci),
                            extra_per_protein = mean(extra), median_rho = median(rho_oof)), by = panel4]
fit <- glm(extra ~ 0 + panel4 + offset(log(obs_loci + 1)), family = quasipoisson, data = per[!is.na(panel4)])
cr[, fisher_p := sapply(panel4, function(pnl) {
  a <- cr[panel4 == pnl]; b <- cr[panel4 != pnl, .(extra_loci = sum(extra_loci), obs_loci = sum(obs_loci))]
  fisher.test(matrix(c(a$extra_loci, a$obs_loci, b$extra_loci, b$obs_loci), 2))$p.value })]
cr[, glm_rate := exp(coef(fit))[paste0("panel4", panel4)]]
setorder(cr, -extra_per_protein); fwrite(cr, paste0(D, "category_rates.csv")); print(cr, digits = 3)

## ---- LaTeX rows of Supplementary Table tab:rho_strata --------------------------
fm <- function(x) format(x, big.mark = ",")
rows <- g[, sprintf("%s & %s & %.3f & %s & %s & %.2f & %s & %s & %.2f & %s / %s \\\\", sub(">=", "$\\\\geq$", sub("<", "$<$", rho_bin)),
                    fm(proteins), median_rho, fm(syn_loci), fm(obs_loci), loci_ratio, fm(syn_snps), fm(obs_snps), snp_ratio,
                    fm(proteins_more_loci), fm(proteins_fewer_loci))]
tot <- per[, sprintf("All & %s & %.3f & %s & %s & %.2f & %s & %s & %.2f & %s / %s \\\\", fm(.N), median(rho_oof), fm(sum(syn_loci)),
                     fm(sum(obs_loci)), sum(syn_loci) / sum(obs_loci), fm(sum(syn_snps)), fm(sum(obs_snps)), sum(syn_snps) / sum(obs_snps),
                     fm(sum(syn_loci > obs_loci)), fm(sum(syn_loci < obs_loci)))]
writeLines(c(rows, "\\midrule", tot), paste0(D, "tab_rho_strata_rows.tex"))

## Reply to Reviewer 1, minor comment (c) / Reviewer 3, minor comment 3 (revision
## 2026-10-20): is the over-representation of the additional Syn-PALM associations among
## cardiometabolic proteins explained by their higher prediction accuracy? Quasi-Poisson
## regression of the additional primary associations per protein on the Olink panel,
## with log(Obs-GWAS primary associations + 1) as offset, without and with a natural
## spline (3 df) of the out-of-fold correlation; likelihood-ratio (F) test of the panel.
##   Rscript code/post_cf3/category_adjusted_rho.R
## Writes results/table_results_cf_noPRS_noPC/category_adjusted_rho.csv
suppressPackageStartupMessages({ library(data.table); library(splines) })
setwd("<PROJECT_DIR>")
D <- "results/table_results_cf_noPRS_noPC/"
per <- fread(paste0(D, "summary_per_protein.csv"))
per[, `:=`(syn = `SynG_loci_5e-08`, obs = `ObsG_loci_5e-08`)]
pm <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
pm[, `:=`(key = gsub("-", "_", id), panel4 = sub("_II$", "", panel))]
per[, key := gsub("-", "_", protein)]   # e.g. hla-e in the results, hla_e in the protein map
per[pm, on = .(key = key), panel4 := i.panel4]
d <- per[!is.na(panel4)][, extra := pmax(syn - obs, 0)]
d[, panel4 := relevel(factor(panel4), ref = "Neurology")]
f0 <- glm(extra ~ offset(log(obs + 1)) + ns(rho_oof, 3), family = quasipoisson, data = d)
f1 <- glm(extra ~ offset(log(obs + 1)) + panel4, family = quasipoisson, data = d)
f2 <- glm(extra ~ offset(log(obs + 1)) + ns(rho_oof, 3) + panel4, family = quasipoisson, data = d)
fnull <- glm(extra ~ offset(log(obs + 1)), family = quasipoisson, data = d)
p_unadj <- anova(fnull, f1, test = "F")$`Pr(>F)`[2]
p_adj   <- anova(f0, f2, test = "F")$`Pr(>F)`[2]
p_rho   <- anova(f1, f2, test = "F")$`Pr(>F)`[2]
rr <- function(f) { co <- summary(f)$coefficients; i <- grep("^panel4", rownames(co))
  data.table(panel = sub("^panel4", "", rownames(co)[i]), rate_ratio = exp(co[i, 1]),
             lo = exp(co[i, 1] - 1.96 * co[i, 2]), hi = exp(co[i, 1] + 1.96 * co[i, 2]), p = co[i, 4]) }
res <- rbind(rr(f1)[, model := "panel only"], rr(f2)[, model := "panel + spline(rho)"])
fwrite(res, paste0(D, "category_adjusted_rho.csv"))
print(res, digits = 3)
cat(sprintf("panel effect: unadjusted P = %.2g; adjusted for rho P = %.2g; rho given panel P = %.2g\n", p_unadj, p_adj, p_rho))
print(d[, .(proteins = .N, median_rho = median(rho_oof), extra_per_protein = mean(extra)), by = panel4])

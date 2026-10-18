## Summaries of the replicate simulations S36 / S4 / S5 (rep_sims_rev1.R).
##   Rscript analyze_rep_rev1.R
## Per study, setting and method, over the causal SNPs (beta != 0):
##   bias, empirical SD of beta-hat, mean model-based SE, power at 5e-8 and 1e-5,
##   and the SE / SD ratios relative to Obs-GWAS (ObsG); over the null SNPs (S5):
##   type I error at 0.05, 1e-3, 1e-4. Writes simulation_code/results_rev1/summary_<study>.csv.
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
R <- "simulation_code/results_rev1/"
for (study in c("S36", "S4", "S5")) {
  files <- list.files(paste0(R, study), "^setting_\\d+_array_\\d+\\.rds$", full.names = TRUE)
  files <- files[as.integer(sub(".*_array_(\\d+)\\.rds", "\\1", files)) <= 500]   # drop test arrays
  if (!length(files)) next
  set_id <- as.integer(sub(".*setting_(\\d+)_.*", "\\1", files))
  tab <- rbindlist(lapply(split(files, set_id), function(fs) {
    x  <- rbindlist(lapply(fs, function(f) { d <- as.data.table(readRDS(f)); d[, file := f] }))
    st <- attr(readRDS(fs[1]), "settings")
    x[, rep_id := .GRP, by = .(file, rep)]
    cz <- x[beta != 0, .(n_est = .N, bias = mean(hat_beta - beta), emp_sd = sd(hat_beta),
                         mean_se = mean(sqrt(var_hat_beta)),
                         power_5e8 = mean(nlp >= -log10(5e-8)), power_1e5 = mean(nlp >= 5),
                         n_lab = mean(n_lab)), by = method]
    nz <- x[beta == 0, .(n_null = .N, t1e_05 = mean(nlp >= -log10(0.05)),
                         t1e_1e3 = mean(nlp >= 3), t1e_1e4 = mean(nlp >= 4)), by = method]
    out <- if (nrow(nz)) merge(cz, nz, by = "method", all = TRUE) else cz
    ref <- out[method == "ObsG"]
    out[, `:=`(se_ratio_vs_ObsG = mean_se / ref$mean_se, sd_ratio_vs_ObsG = emp_sd / ref$emp_sd,
               power_ratio_vs_ObsG = power_5e8 / ref$power_5e8)]
    keep <- setdiff(names(st), c("Gm", "Em"))
    cbind(setting = as.integer(sub(".*setting_(\\d+)_.*", "\\1", fs[1])), n_files = length(fs),
          as.data.table(lapply(st[keep], function(v) if (length(v) == 1) v else paste(v, collapse = ";"))), out)
  }), fill = TRUE)
  fwrite(tab, sprintf("%ssummary_%s.csv", R, study))
  cat("\n==", study, ":", uniqueN(tab$setting), "settings\n")
  print(tab[method %in% c("SynSurrG", "ObsG", "SynSurr"),
            intersect(c("setting", "rho_g", "rho_e", "miss", "prop_related", "label", "method", "n_est", "bias",
                        "emp_sd", "mean_se", "se_ratio_vs_ObsG", "power_5e8", "t1e_1e3"), names(tab)), with = FALSE],
        digits = 4)
}

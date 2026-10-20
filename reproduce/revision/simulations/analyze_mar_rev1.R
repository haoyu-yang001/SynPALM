## Summaries of rep_sims_mar_rev1.R (S5W: labelling on W; S5V: S5 data refitted with true VCs).
##   Rscript analyze_mar_rev1.R
## Per setting and method, over the causal SNPs: bias, relative bias (bias / beta_G) with its
## Monte Carlo SE (from the replicate means, which allows for the dependence of the 20 causal
## SNPs of a replicate), empirical SD, mean SE, power at 5e-8; over the null SNPs: type I error
## at 1e-3. Also the mean method-of-moments VCs and the realised corr(W, Y | Yhat).
## Checks: S5V must reproduce results_rev1/S5/ exactly for the MoM methods; MCAR rows of S5W
## are compared with S5 MCAR within Monte Carlo error.
## Writes simulation_code/results_rev1/S5/summary_mar_w.csv, summary_true_vc.csv.
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
R <- "simulation_code/results_rev1/S5/"
nrep_expect <- 25L; narr_expect <- 20L

read_set <- function(fs) {
  x <- rbindlist(lapply(fs, function(f) { d <- as.data.table(readRDS(f)); d[, file := basename(f)] }), fill = TRUE)
  x[, rep_id := .GRP, by = .(file, rep)]; x
}
summ <- function(x, beta_G) {
  rb <- x[beta != 0, .(m = mean(hat_beta - beta)), by = .(method, rep_id)][
          , .(mc_se_bias = sd(m) / sqrt(.N), n_rep = .N), by = method]
  cz <- x[beta != 0, .(n_est = .N, bias = mean(hat_beta - beta), emp_sd = sd(hat_beta),
                       mean_se = mean(sqrt(var_hat_beta)), power_5e8 = mean(nlp >= -log10(5e-8)),
                       n_lab = mean(n_lab), tau_T2 = mean(tau_T2), tau_TS = mean(tau_TS), tau_S2 = mean(tau_S2),
                       sigma_T2 = mean(sigma_T2), sigma_TS = mean(sigma_TS), sigma_S2 = mean(sigma_S2),
                       pcor_W = mean(pcor_W_Y_given_Yhat), pcor_W_x = mean(pcor_W_Y_given_Yhat_x)), by = method]
  nz <- x[beta == 0, .(n_null = .N, t1e_1e3 = mean(nlp >= 3)), by = method]
  out <- merge(merge(cz, rb, by = "method"), nz, by = "method")
  out[, `:=`(rel_bias = bias / beta_G, rel_bias_mcse = mc_se_bias / beta_G)]
  out[, se_ratio_vs_ObsG := mean_se / out[method == "ObsG", mean_se]]
  out[]
}

res <- list()
for (sub in c("mar_w", "true_vc")) {
  files <- list.files(paste0(R, sub), "^setting_\\d+_array_\\d+\\.rds$", full.names = TRUE)
  files <- files[as.integer(sub(".*_array_(\\d+)\\.rds", "\\1", files)) <= 500]   # drop test arrays
  if (!length(files)) next
  set_id <- as.integer(sub(".*setting_(\\d+)_.*", "\\1", files))
  tab <- rbindlist(lapply(split(files, set_id), function(fs) {
    x <- read_set(fs); st <- attr(readRDS(fs[1]), "settings")
    if (length(fs) != narr_expect || x[, uniqueN(rep_id)] != narr_expect * nrep_expect)
      warning(sprintf("%s setting %s: %d files, %d replicates", sub, st$label, length(fs), x[, uniqueN(rep_id)]))
    s <- summ(x, st$beta_G)
    if (sub == "true_vc") {
      ## exact reproduction of results_rev1/S5/ and paired difference trueVC - MoM
      orig <- read_set(file.path(R, basename(fs)))
      key <- c("file", "rep", "method", "snp")
      mm <- merge(x[, c(key, "hat_beta", "nlp"), with = FALSE], orig[, c(key, "hat_beta", "nlp"), with = FALSE],
                  by = key, suffixes = c("", "_S5"))
      s <- merge(s, mm[, .(max_abs_diff_beta_vs_S5 = max(abs(hat_beta - hat_beta_S5)),
                           max_abs_diff_nlp_vs_S5 = max(abs(nlp - nlp_S5)), n_matched = .N), by = method],
                 by = "method", all.x = TRUE)
      pd <- dcast(x[beta != 0, .(m = mean(hat_beta - beta)), by = .(method, rep_id)], rep_id ~ method, value.var = "m")
      s[method == "SynSurrG_trueVC", paired_diff_vs_MoM := mean(pd$SynSurrG_trueVC - pd$SynSurrG) / st$beta_G]
      s[method == "SynSurrG_trueVC", paired_diff_mcse := sd(pd$SynSurrG_trueVC - pd$SynSurrG) / sqrt(nrow(pd)) / st$beta_G]
      s[method == "SynSurr_trueVC", paired_diff_vs_MoM := mean(pd$SynSurr_trueVC - pd$SynSurr) / st$beta_G]
      s[method == "SynSurr_trueVC", paired_diff_mcse := sd(pd$SynSurr_trueVC - pd$SynSurr) / sqrt(nrow(pd)) / st$beta_G]
    }
    cbind(setting = as.integer(sub(".*setting_(\\d+)_.*", "\\1", fs[1])), n_files = length(fs),
          label = st$label, r_w = if (is.null(st$r_w)) NA_real_ else st$r_w, miss = st$miss,
          beta_G = st$beta_G, seed_first = st$seed - st$array + 1L, s)
  }), fill = TRUE)
  fwrite(tab, sprintf("%ssummary_%s.csv", R, sub))
  res[[sub]] <- tab
  cat("\n==", sub, "\n")
  print(tab[, intersect(c("setting", "label", "r_w", "miss", "method", "n_rep", "rel_bias", "rel_bias_mcse",
                          "power_5e8", "t1e_1e3", "se_ratio_vs_ObsG", "pcor_W", "pcor_W_x", "max_abs_diff_beta_vs_S5",
                          "paired_diff_vs_MoM", "paired_diff_mcse"), names(tab)), with = FALSE], digits = 3)
  print(unique(tab[, .(setting, label, r_w, miss, n_lab, tau_T2, tau_TS, tau_S2, sigma_T2, sigma_TS, sigma_S2)]), digits = 3)
}

## MCAR rows of S5W vs Table S5 (results_rev1/summary_S5.csv), difference in units of its MC SE
if (!is.null(res$mar_w) && file.exists("simulation_code/results_rev1/summary_S5.csv")) {
  s5 <- fread("simulation_code/results_rev1/summary_S5.csv")[label == "random"]
  w  <- res$mar_w[label == "random"]
  cmp <- merge(w[, .(miss, method, bias_new = bias, sd_new = emp_sd, n_new = n_est, pow_new = power_5e8, t1e_new = t1e_1e3)],
               s5[, .(miss, method, bias_S5 = bias, sd_S5 = emp_sd, n_S5 = n_est, pow_S5 = power_5e8, t1e_S5 = t1e_1e3)],
               by = c("miss", "method"))
  cmp[, `:=`(z_bias = (bias_new - bias_S5) / sqrt(sd_new^2 / n_new + sd_S5^2 / n_S5),
             z_power = (pow_new - pow_S5) / sqrt(pow_new * (1 - pow_new) / n_new + pow_S5 * (1 - pow_S5) / n_S5),
             z_t1e = (t1e_new - t1e_S5) / sqrt(1e-3 * (1 - 1e-3) * 2 / 1e5))]
  cat("\n== MCAR check: S5W random vs S5 random\n")
  print(cmp[, .(miss, method, bias_new, bias_S5, z_bias, pow_new, pow_S5, z_power, t1e_new, t1e_S5, z_t1e)], digits = 3)
}

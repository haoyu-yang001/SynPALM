## Inputs for the power (Fig 3, S5-S7) and bias (Fig 2A) figures of the revision, in the
## file names and formats the submitted plotting scripts read. All rho are rerun with
## power_bias_rho075_rev1.R (rho = 0.75 first; rho = 0, 0.25, 0.5 from 2026-10-16); a
## setting whose rerun is complete is merged from the rerun arrays, otherwise the
## submitted result is used through a symbolic link (the link is removed, never written
## through, when the merged rerun file replaces it):
##   simulation_code/results_rev1/power_kinsetting_results/  (6 columns of -log10 p)
##   simulation_code/results_rev1/point_estimate/           (18 columns)
##   Rscript merge_power_bias_rev1.R
setwd("<PROJECT_DIR>")
R0 <- "simulation_code/results/"; R1 <- "simulation_code/results_rev1/"
for (d in c("power_kinsetting_results", "point_estimate")) {
  out <- paste0(R1, d, "/"); dir.create(out, showWarnings = FALSE)
  old <- list.files(paste0(R0, d), "\\.rds$")
  old <- old[!grepl("_rho_0.75_", old)]
  for (f in old) if (!file.exists(paste0(out, f))) file.symlink(normalizePath(paste0(R0, d, "/", f)), paste0(out, f))
}
kin <- c("0.5" = "siblings", "0.25" = "halfsiblings", "0.125" = "cousin")
miss_vec <- c(0, 0.25, 0.5, 0.75, 0.9)
report <- NULL
for (mode in c("power", "bias")) for (rho in c(0, 0.25, 0.5, 0.75)) {
  if (mode == "bias" && !rho %in% c(0, 0.75)) next
  src <- sprintf("%s%s_rho%03d/", R1, mode, round(100 * rho))
  grid <- if (mode == "power") expand.grid(GRMr = c(0.125, 0.25, 0.5), miss = miss_vec, her = seq(0.0001, 0.001, 0.0001))
          else expand.grid(GRMr = 0.5, miss = miss_vec, her = 0)
  n_arr <- if (mode == "power") 200 else 1000
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]; k <- kin[[as.character(g$GRMr)]]
    fs <- sprintf("%sGRMr_%s_missing_%s_rho_%s_rhogeno_%s_her_%s_chunkarray_%d_rev1.rds", src, g$GRMr, g$miss, rho, k, g$her, seq_len(n_arr))
    have <- file.exists(fs)
    report <- rbind(report, data.frame(mode, rho, GRMr = g$GRMr, miss = g$miss, her = g$her, arrays = sum(have)))
    if (!all(have)) next
    x <- do.call(rbind, lapply(fs, readRDS))
    if (mode == "power") {
      x <- x[, paste0(c("OracleG", "ObsG", "SynSurrG", "Oracle", "Obs", "SynSurr"), "_negative_log10_pval")]
      colnames(x) <- sub("_negative_log10_pval", "", colnames(x))
      f_out <- sprintf("%spower_kinsetting_results/GRMr_%s_missing_%s_rho_%s_rhogeno_%s_her_%s_negative_log10_pvals_for_power.rds",
                       R1, g$GRMr, g$miss, rho, k, g$her)
      if (nzchar(Sys.readlink(f_out))) file.remove(f_out)      # drop the link to the submitted file
      saveRDS(x, f_out)
    } else {
      f_out <- sprintf("%spoint_estimate/GRMr_%s_missing_%s_rho_%s_rhogeno_%s_her_0_negative_log10_pvals_for_power.rds",
                       R1, g$GRMr, g$miss, rho, k)
      if (nzchar(Sys.readlink(f_out))) file.remove(f_out)
      saveRDS(x, f_out)
    }
  }
}
print(aggregate(arrays ~ mode + rho, report, function(a) c(complete = sum(a == max(a)), min = min(a))))
cat("settings complete:", sum(report$arrays == ifelse(report$mode == "power", 200, 1000)), "of", nrow(report), "\n")

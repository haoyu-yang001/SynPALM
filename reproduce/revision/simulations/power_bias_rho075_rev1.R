## Rerun of the rho = 0.75 settings of the submitted power (Fig 3, Supp Figs
## S5-S7; simulation_power_v2.R) and bias (Fig 2A;
## simulation_point_estimate_for_robustness.R) simulations, whose genetic
## covariance broke Cauchy-Schwarz (Reviewer 2, comment 5). Same data generation,
## methods and output columns as the submitted scripts; changes:
##   - rho = 0.75 uses sim_cov_TS() (genetic correlation 0.99, total rho 0.75);
##   - every replicate block is seeded;
##   - the genotype model matches the GRM (r = 0.5 siblings, 0.25 half-siblings,
##     0.125 cousins), as in the submitted power runs.
##   Rscript power_bias_rho075_rev1.R power <GRMr> <array> [n_snps=50] [rho=0.75]   (200 arrays x 50 SNPs)
##   Rscript power_bias_rho075_rev1.R bias 0.5 <array> [n_snps=50] [rho=0.75]       (1000 arrays x 50 SNPs)
## Revision 2 (2026-10-16): the settings rho = 0, 0.25, 0.5 kept from the submission are
## also rerun with this script (argument rho), so that all power and bias results come
## from the released code with fixed seeds; output in results_rev1/<mode>_rho<100 rho>/,
## seeds offset by 1e7 x (1, 2, 3) for rho = 0, 0.25, 0.5 (rho = 0.75 unchanged).
setwd("<PROJECT_DIR>")
suppressMessages(source("code/funcs_v2.R"))
options(synsurrg.verbose = FALSE)
library(withr)
args     <- commandArgs(trailingOnly = TRUE)
mode     <- args[1]
GRMr     <- as.numeric(args[2])
array_id <- as.integer(args[3])
kinsetting <- c("0.5" = "siblings", "0.25" = "halfsiblings", "0.125" = "cousin")[[as.character(GRMr)]]
num_snps_each_array <- if (length(args) >= 4) as.integer(args[4]) else 50L
n_obs <- 3e4
miss_vec <- c(0, 0.25, 0.5, 0.75, 0.9)
rho <- if (length(args) >= 5) as.numeric(args[5]) else 0.75
rho_off <- c("0.75" = 0L, "0" = 10000000L, "0.25" = 20000000L, "0.5" = 30000000L)[[as.character(rho)]]
her_vec <- if (mode == "power") seq(0.0001, 0.001, 0.0001) else 0
tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5; pve_x <- 0; maf <- 0.25
OUT <- sprintf("simulation_code/results_rev1/%s_rho%03d/", mode, round(100 * rho))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cols <- c("OracleG_negative_log10_pval","OracleG_hat_beta","OracleG_var_hat_beta",
          "ObsG_negative_log10_pval","ObsG_hat_beta","ObsG_var_hat_beta",
          "SynSurrG_negative_log10_pval","SynSurrG_hat_beta","SynSurrG_var_hat_beta",
          "Oracle_negative_log10_pval","Oracle_hat_beta","Oracle_var_hat_beta",
          "Obs_negative_log10_pval","Obs_hat_beta","Obs_var_hat_beta",
          "SynSurr_negative_log10_pval","SynSurr_hat_beta","SynSurr_var_hat_beta")

for (mi in seq_along(miss_vec)) {
  miss <- miss_vec[mi]
  GRM <- DGP_GRMr(n_obs, miss, GRMr)
  hh <- find_blocks_vectorized(GRM)
  independent_indices <- with_seed(25324, sapply(hh, function(x) x[sample.int(length(x), 1)]))
  pre1 <- DGP_GRMr_nullmodel_step1(n_obs, miss, tauT, tauS, sigmaT, sigmaS, GRM)
  pre2 <- DGP_GRMr_nullmodel_step2(n_obs, miss, tauT, tauS, sigmaT, sigmaS, rho, GRM, pre1)
  for (hi in seq_along(her_vec)) {
    pve_g <- her_vec[hi]
    f <- sprintf("%sGRMr_%s_missing_%s_rho_%s_rhogeno_%s_her_%s_chunkarray_%d_rev1.rds",
                 OUT, GRMr, miss, rho, kinsetting, pve_g, array_id)
    if (file.exists(f)) next
    seed <- 20261015L + rho_off + (if (mode == "power") 0L else 5000000L) + 1000000L * match(GRMr, c(0.125, 0.25, 0.5)) +
            1000L * array_id + 100L * mi + hi
    set.seed(seed)
    cat("miss", miss, "rho", rho, "pve_g", pve_g, "array", array_id, "GRMr", GRMr, "seed", seed, "\n")
    res <- matrix(NA_real_, num_snps_each_array, length(cols), dimnames = list(NULL, cols))
    for (ss in seq_len(num_snps_each_array)) {
      mydf <- DGP_GRMr_gmodel_step3(n_obs, miss, tauT = tauT, sigmaT = sigmaT, pve_x = pve_x, pve_g = pve_g,
                                    maf, kinsetting, Genotype_each = NULL, pre2, GRM)
      G_all <- mydf$G_all
      r1 <- score_test_OracleG_single(G_all, step1_pars = OracleG_typeIerror_step1(mydf))
      r2 <- score_test_ObsG_single(G_all, step1_pars = ObsG_typeIerror_step1(mydf))
      r3 <- score_test_SynSurrG_single(G_all, step1_pars = SynSurrG_typeIerror_step1(mydf))
      r4 <- score_test_Oracle_single(G_all, step1_pars = Oracle_typeIerror_step1(mydf, independent_indices), independent_indices)
      r5 <- score_test_Obs_single(G_all, step1_pars = Obs_typeIerror_step1(mydf, independent_indices), independent_indices)
      r6 <- score_test_SynSurr_single(G_all, step1_pars = SynSurr_typeIerror_step1(mydf, independent_indices), independent_indices)
      res[ss, ] <- c(r1$negative_log10_pval_OracleG, r1$hat_beta_OracleG, r1$var_hat_beta_OracleG,
                     r2$negative_log10_pval_ObsG, r2$hat_beta_ObsG, r2$var_hat_beta_ObsG,
                     r3$negative_log10_pval_SynSurrG, r3$hat_beta_SynSurrG, r3$var_hat_beta_SynSurrG,
                     r4$negative_log10_pval_Oracle, r4$hat_beta_Oracle, r4$var_hat_beta_Oracle,
                     r5$negative_log10_pval_Obs, r5$hat_beta_Obs, r5$var_hat_beta_Obs,
                     r6$negative_log10_pval_SynSurr, r6$hat_beta_SynSurr, r6$var_hat_beta_SynSurr)
    }
    attr(res, "settings") <- list(mode = mode, GRMr = GRMr, kinsetting = kinsetting, miss = miss, rho = rho,
                                  pve_g = pve_g, seed = seed, cov_TS = sim_cov_TS(rho))
    saveRDS(res, f)
  }
}

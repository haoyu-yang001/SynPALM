## Type I error when relatives are split between the labelled and the unlabelled
## sets (revision; Reviewer comments on POP-GWAS and relatedness). Same model,
## parameters and tests as typeIerror_rev1.R for sib pairs (r = 0.5, n_obs =
## 30,000 labelled, tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5, MAF 0.25,
## no covariate effect, rho = corr(Y, Yhat) with the 3/5 : 2/5 split of
## sim_cov_TS), except that the labelled individuals are drawn completely at
## random (label_set(type = "random") of sim_dgp_rev1.R), so that for a share of
## the pairs one sibling is labelled and the other is not. POP-GWAS assumes that
## the labelled and the unlabelled samples are independent, which this design
## violates; Syn-PALM models the cross-set relatedness through the GRM.
##   Rscript typeIerror_split_rev1.R <array> [n_chunks=250] [miss list, e.g. 0.5:0.75]
## 250 chunks x 200 SNPs = 5e4 null SNPs per array and (pi, rho); 200 arrays give
## 1e7 SNPs per setting. Output in the format of typeIerror_rev1.R, so that
## analyze_t1e_rev1.R (T1E_DIR) summarises it.
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R"); source("pop_gwas/Fun.R")
                   source("simulation_code/rev1/sim_dgp_rev1.R") })
options(synsurrg.verbose = FALSE)
args       <- commandArgs(trailingOnly = TRUE)
array_id   <- as.integer(args[1])
num_chunks <- if (length(args) >= 2) as.integer(args[2]) else 250L
miss_vec   <- if (length(args) >= 3) as.numeric(strsplit(args[3], "[,:]")[[1]]) else c(0.5, 0.75)
rho_vec    <- c(0, 0.25, 0.5, 0.75)
n_obs <- 3e4; maf <- 0.25; chunk_size <- 200L
OUT <- "simulation_code/results_rev1/typeIerror_split/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

for (miss in miss_vec) {
  out_file <- sprintf("%sGRMr_0.5_split_missing_%s_array_%d.rds", OUT, miss, array_id)
  if (file.exists(out_file)) { tick("exists: ", out_file); next }
  seed <- 20261101L + 1000L * match(miss, c(0.25, 0.5, 0.75, 0.9)) + array_id
  set.seed(seed)
  n_all <- round(n_obs / (1 - miss)); n_miss <- n_all - n_obs
  coh <- make_cohort(n_all, prop_related = 1, r = 0.5)
  lab <- label_set(coh, miss, "random")
  x <- rnorm(n_all)
  split_pairs <- mean((coh$pair1 %in% lab) != (coh$pair2 %in% lab))
  fits <- lapply(rho_vec, function(rho) {
    e <- draw_noise(coh, vc_from_total(rho))
    mydf <- make_mydf(coh, Y = e[, 1], S = e[, 2], x = x, lab = lab)
    list(syn = add_fastscan_cache(SynSurrG_typeIerror_step1(mydf)), obs = ObsG_typeIerror_step1(mydf),
         obs_hat = ObsG_hatY_typeIerror_step1(mydf), unobs_hat = unObsG_hatY_typeIerror_step1(mydf),
         obs_lm = ObsG_typeIerror_step1_lm(mydf), unobs_lm = unObsG_hatY_typeIerror_step1_lm(mydf))
  })
  tick("miss ", miss, ": null models fitted (seed ", seed, "; ", round(100 * split_pairs, 1), "% of pairs split)")
  res <- vector("list", length(rho_vec))
  for (ci in seq_len(num_chunks)) {
    G <- gene_drop_cohort(coh, chunk_size, maf)
    for (ri in seq_along(rho_vec)) {
      f <- fits[[ri]]
      s  <- score_test_SynSurrG_fast(G, f$syn)
      o  <- score_test_ObsG_multiply(g_matrix = G, step1_pars = f$obs)
      h1 <- score_test_ObsG_hatY_multiply(g_matrix = G, step1_pars = f$obs_hat)
      h2 <- score_test_unObsG_hatY_multiply(g_matrix = G, step1_pars = f$unobs_hat)
      l1 <- score_test_ObsG_lm_multiply(g_matrix = G, step1_pars = f$obs_lm)
      l2 <- score_test_unObsG_hatY_lm_multiply(g_matrix = G, step1_pars = f$unobs_lm)
      pop <- cal_qt_opt_vec(
        z_df = data.frame(Z1 = h1$hat_beta_ObsG_hatY / sqrt(h1$var_hat_beta_ObsG_hatY),
                          Z2 = h2$hat_beta_unObsG_hatY / sqrt(h2$var_hat_beta_unObsG_hatY),
                          Z3 = o$hat_beta_ObsG / sqrt(o$var_hat_beta_ObsG)),
        n = n_obs * l1$var_hat_beta_ObsG / o$var_hat_beta_ObsG,
        N = n_miss * l2$var_hat_beta_unObsG_hatY / h2$var_hat_beta_unObsG_hatY,
        r = rho_vec[ri], eaf = maf)
      res[[ri]][[ci]] <- cbind(SynSurrG_negative_log10_pvalues = s$negative_log10_pval_SynSurrG,
                               ObsG_negative_log10_pvalues = o$negative_log10_pval_ObsG, POPGWAS_p = pop$P)
    }
    if (ci %% 50 == 0) tick("chunk ", ci, "/", num_chunks)
  }
  out <- setNames(lapply(res, function(x) do.call(rbind, x)), paste0("rho_", rho_vec))
  attr(out, "settings") <- list(GRMr = 0.5, design = "random labelling (split pairs)", miss = miss, n_obs = n_obs,
                                n_all = n_all, split_pairs = split_pairs, seed = seed, num_chunks = num_chunks, maf = maf)
  saveRDS(out, out_file)
  tick("saved ", out_file)
}

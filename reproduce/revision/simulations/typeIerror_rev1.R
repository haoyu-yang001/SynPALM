## Type I error simulations for the revision (Table 1, Fig 2B, Supp Figs S2-S4).
## Replaces simulation_code/simluation_typeIerror_control_v2.R and
## pop_gwas/simulation_error_control_{1..4}.R:
##   - every random draw is seeded (seed printed and saved with the output);
##   - the genotype model matches the GRM: r = 0.5 siblings, 0.25 half-siblings,
##     0.125 first cousins (the submitted runs used half-siblings throughout);
##   - rho = 0.75 uses the corrected covariances (sim_cov_TS in funcs_v2.R);
##   - Syn-PALM, Obs-GWAS and (sib pairs only) POP-GWAS are scored on the same
##     null phenotype and the same SNPs; POP-GWAS at pi = 0 has no unlabelled
##     individuals and reduces to the labelled-only LMM test (Obs-GWAS);
##   - the same simulated SNPs are scored against the four rho phenotypes of a
##     missing rate (one null phenotype per (pi, rho) per array, as before);
##   - Syn-PALM is scored with the vectorised test of funcs_fastscan.R (identical
##     to score_test_SynSurrG_multiply up to 1e-13), and the gene drop is done for
##     a whole chunk of SNPs at once (same model as generate_gmatrix_kinship);
##   - Obs-GWAS and POP-GWAS are computed for sib pairs (r = 0.5) only.
## Data generation (unchanged, funcs_v2.R DGP_GRMr_*): block-diagonal GRM of
## pairs with relatedness r; S ~ N(0, tauS^2 K + sigmaS^2 I); Y_obs | S from the
## conditional normal with Cov(Y, S) = tauTS K + sigmaTS I; the first n_obs =
## 30,000 individuals are labelled; SNPs by gene drop with MAF 0.25, no effect.
##   Rscript typeIerror_rev1.R <GRMr> <array> [n_chunks=1000] [miss list, e.g. 0,0.5,0.75]
## 1000 chunks x 200 SNPs = 2e5 null SNPs per array. Sib pairs: 500 arrays
## (1e8 SNPs per (pi, rho)) for pi = 0, 0.5, 0.75 (Table 1), 100 arrays (2e7)
## for pi = 0.25, 0.9 (QQ plots only); half-siblings and cousins: 100 arrays.
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R"); source("pop_gwas/Fun.R") })
options(synsurrg.verbose = FALSE)

args       <- commandArgs(trailingOnly = TRUE)
GRMr       <- as.numeric(args[1])
array_id   <- as.integer(args[2])
num_chunks <- if (length(args) >= 3) as.integer(args[3]) else 1000L
kinsetting <- c("0.5" = "siblings", "0.25" = "halfsiblings", "0.125" = "cousin")[[as.character(GRMr)]]
gi         <- match(GRMr, c(0.125, 0.25, 0.5))
do_pop     <- GRMr == 0.5

miss_all <- c(0, 0.25, 0.5, 0.75, 0.9)
miss_vec <- if (length(args) >= 4) as.numeric(strsplit(args[4], ",")[[1]]) else miss_all
rho_vec  <- c(0, 0.25, 0.5, 0.75)
n_obs <- 3e4; tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5
pve_x <- 0; maf <- 0.25; chunk_size <- 200

## gene drop for m SNPs at once: pair (2i-1, 2i) shares grandparents
## (siblings: children of the same two parents; half-siblings: one shared
## parent; cousins: children of two siblings), as in generate_gmatrix_kinship()
gene_drop <- function(n_all, m, maf, kinsetting) {
  np <- n_all / 2
  founder <- function() matrix(rbinom(np * m, 2L, maf), np, m)
  gam <- function(g) matrix(rbinom(length(g), 1L, g / 2), nrow(g))
  gf <- founder(); gm <- founder()
  pair <- switch(kinsetting,
    siblings     = list(gam(gf) + gam(gm), gam(gf) + gam(gm)),
    halfsiblings = list(gam(gf) + gam(gm), gam(gf) + gam(founder())),
    cousin       = { f <- gam(gf) + gam(gm); u <- gam(gf) + gam(gm)
                     list(gam(f) + gam(founder()), gam(u) + gam(founder())) })
  G <- matrix(0L, n_all, m)
  G[seq(1, n_all, 2), ] <- pair[[1]]; G[seq(2, n_all, 2), ] <- pair[[2]]
  G
}

OUT <- "simulation_code/results_rev1/typeIerror/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

for (miss in miss_vec) {
  mi <- match(miss, miss_all)
  out_file <- sprintf("%sGRMr_%s_%s_missing_%s_array_%d.rds", OUT, GRMr, kinsetting, miss, array_id)
  if (file.exists(out_file)) { tick("exists: ", out_file); next }
  seed <- 20261013L + 100000L * gi + 100L * array_id + mi
  set.seed(seed)
  n_miss <- round(n_obs * miss / (1 - miss))
  GRM <- DGP_GRMr(n_obs, miss, GRMr)
  pre1 <- DGP_GRMr_nullmodel_step1(n_obs, miss, tauT, tauS, sigmaT, sigmaS, GRM)

  ## one null phenotype per rho, null models fitted once
  fits <- lapply(rho_vec, function(rho) {
    pre2 <- DGP_GRMr_nullmodel_step2(n_obs, miss, tauT, tauS, sigmaT, sigmaS, rho, GRM, pre1)
    mydf <- DGP_GRMr_nullmodel_step3(n_obs, miss, tauT, sigmaT, pve_x, pre2, GRM)
    f <- list(syn = add_fastscan_cache(SynSurrG_typeIerror_step1(mydf)))
    if (do_pop) f$obs <- ObsG_typeIerror_step1(mydf)
    if (do_pop && miss > 0) {
      f$obs_hat   <- ObsG_hatY_typeIerror_step1(mydf)
      f$unobs_hat <- unObsG_hatY_typeIerror_step1(mydf)
      f$obs_lm    <- ObsG_typeIerror_step1_lm(mydf)
      f$unobs_lm  <- unObsG_hatY_typeIerror_step1_lm(mydf)
    }
    f
  })
  rm(pre1); invisible(gc())
  tick("GRMr ", GRMr, " miss ", miss, ": null models fitted (seed ", seed, ")")

  res <- vector("list", length(rho_vec))
  for (ci in seq_len(num_chunks)) {
    G <- gene_drop(n_obs + n_miss, chunk_size, maf, kinsetting)
    for (ri in seq_along(rho_vec)) {
      f <- fits[[ri]]
      s <- score_test_SynSurrG_fast(G, f$syn)
      blk <- cbind(SynSurrG_negative_log10_pvalues = s$negative_log10_pval_SynSurrG)
      if (do_pop) {
        o <- score_test_ObsG_multiply(g_matrix = G, step1_pars = f$obs)
        blk <- cbind(blk, ObsG_negative_log10_pvalues = o$negative_log10_pval_ObsG)
        if (miss > 0) {
          z3 <- o$hat_beta_ObsG / sqrt(o$var_hat_beta_ObsG)
          h1 <- score_test_ObsG_hatY_multiply(g_matrix = G, step1_pars = f$obs_hat)
          h2 <- score_test_unObsG_hatY_multiply(g_matrix = G, step1_pars = f$unobs_hat)
          l1 <- score_test_ObsG_lm_multiply(g_matrix = G, step1_pars = f$obs_lm)
          l2 <- score_test_unObsG_hatY_lm_multiply(g_matrix = G, step1_pars = f$unobs_lm)
          pop <- cal_qt_opt_vec(
            z_df = data.frame(Z1 = h1$hat_beta_ObsG_hatY / sqrt(h1$var_hat_beta_ObsG_hatY),
                              Z2 = h2$hat_beta_unObsG_hatY / sqrt(h2$var_hat_beta_unObsG_hatY),
                              Z3 = z3),
            n = n_obs * l1$var_hat_beta_ObsG / o$var_hat_beta_ObsG,
            N = n_miss * l2$var_hat_beta_unObsG_hatY / h2$var_hat_beta_unObsG_hatY,
            r = rho_vec[ri], eaf = maf)
          blk <- cbind(blk, POPGWAS_p = pop$P)
        } else {
          blk <- cbind(blk, POPGWAS_p = 10^(-o$negative_log10_pval_ObsG))
        }
      }
      res[[ri]][[ci]] <- blk
    }
    if (ci %% 100 == 0) tick("chunk ", ci, "/", num_chunks)
  }
  out <- setNames(lapply(res, function(x) do.call(rbind, x)), paste0("rho_", rho_vec))
  attr(out, "settings") <- list(GRMr = GRMr, kinsetting = kinsetting, miss = miss, rho = rho_vec,
                                n_obs = n_obs, n_miss = n_miss, seed = seed, array = array_id,
                                n_snps = num_chunks * chunk_size,
                                cov_TS = t(sapply(rho_vec, sim_cov_TS)),
                                vc = lapply(fits, function(f) f$syn$vc))
  saveRDS(out, out_file)
  tick("saved ", out_file)
  rm(fits, res, out, GRM); invisible(gc())
}

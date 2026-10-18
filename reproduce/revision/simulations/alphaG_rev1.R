## S2: SNPs that affect the synthetic phenotype but not the target phenotype
## (alpha_G != 0, beta_G = 0; Reviewer 2 comments 2, 4 and 8).
## Sib pairs, n_obs = 30,000 labelled (whole families), rho = 0.5 (3/5 genetic,
## 2/5 residual). Each phenotype pair carries K "alpha SNPs", each explaining h_a
## of the variance of S and nothing of Y, plus 1,000 null SNPs:
##   Y = x b_x + u_T + e_T,   S = sum_j alpha g_j + x b_x + u_S + e_S.
## Every SNP is tested for H0: beta_G = 0 with
##   SynSurrG  Syn-PALM, variance components estimated once under the global null
##             (as in the paper)
##   SynSurrG_refit  Syn-PALM with the variance components re-estimated with the
##             tested SNP as a covariate (first n_refit alpha SNPs only)
##   ObsG, SynSurr  for reference.
## With the variance components fixed, the Syn-PALM statistic of SNP g does not
## change when alpha g is added to S (alpha_G is estimated for every SNP), so any
## inflation could only come through the variance components; the script also
## records max |stat(S) - stat(S - alpha g)| for the first SNPs as a check.
##   Rscript alphaG_rev1.R <setting> <array> [n_rep=20]
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R")
                   source("simulation_code/rev1/sim_dgp_rev1.R") })
options(synsurrg.verbose = FALSE)
args     <- commandArgs(trailingOnly = TRUE)
setting  <- as.integer(args[1])
array_id <- as.integer(args[2])
n_rep    <- if (length(args) >= 3) as.integer(args[3]) else 20L

S <- data.frame(miss = c(0.5, 0.5, 0.5, 0.5, 0.9, 0.9),
                h_a  = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-4, 1e-2),
                K    = c(10000L, 1000L, 100L, 10L, 1000L, 10L))
st <- S[setting, ]
maf <- 0.25; n_obs0 <- 3e4; n_null <- 1000L; n_refit <- 10L; pve_x <- 0.10
vc <- vc_from_total(0.5)
n_all <- round(n_obs0 / (1 - st$miss))
coh <- make_cohort(n_all, prop_related = 1, r = 0.5)
V_S <- vc$Gm[2, 2] + vc$Em[2, 2]; V_T <- vc$Gm[1, 1] + vc$Em[1, 1]
alpha <- sqrt(st$h_a * V_S / (2 * maf * (1 - maf)))

OUT <- "simulation_code/results_rev1/S2/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("%ssetting_%02d_array_%03d.rds", OUT, setting, array_id)
if (file.exists(out_file)) quit(save = "no")
seed <- 20261017L + 10000L * setting + array_id
set.seed(seed)
cat(sprintf("S2 setting %d: miss %s, h_a %g (alpha %.4f), K %d | n %d | seed %d\n",
            setting, st$miss, st$h_a, alpha, st$K, n_all, seed))

## Syn-PALM null model with given (already constrained) variance components
txt <- deparse(SynSurrG_typeIerror_step1)
txt[1] <- sub("function (mydf) ", "function (mydf, vc_fixed = NULL) ", txt[1], fixed = TRUE)
i <- grep("vc <- .vc$vc", txt, fixed = TRUE)[1]
txt <- append(txt, c("if (!is.null(vc_fixed)) {",
  "hat_tau_T2 <- vc_fixed[['tau_T2']]; hat_tau_TS <- vc_fixed[['tau_TS']]; hat_tau_S2 <- vc_fixed[['tau_S2']]",
  "hat_sigma_T2 <- vc_fixed[['sigma_T2']]; hat_sigma_TS <- vc_fixed[['sigma_TS']]; hat_sigma_S2 <- vc_fixed[['sigma_S2']]",
  "vc$used <- unlist(vc_fixed[c('tau_T2','tau_TS','tau_S2','sigma_T2','sigma_TS','sigma_S2')]) }"), after = i)
step1_fixed <- eval(parse(text = txt))
stopifnot(length(i) == 1)

## method-of-moments variance components with covariates W (n x q incl. intercept)
GRM <- coh$GRM
mom_vc <- function(Y_obs, S, W) {
  lab <- which(!is.na(Y_obs))
  rT <- lm.fit(W[lab, , drop = FALSE], Y_obs[lab])$residuals; rS <- lm.fit(W, S)$residuals
  Koo <- GRM[lab, lab]; Go <- Koo; diag(Go) <- 0; Ga <- GRM; diag(Ga) <- 0
  tT <- as.numeric(t(rT) %*% Go %*% rT) / sum(Go^2); tS <- as.numeric(t(rS) %*% Ga %*% rS) / sum(Ga^2)
  tX <- as.numeric(t(rS[lab]) %*% Go %*% rT) / sum(Go^2)
  constrain_vc_bivariate(tT, tX, tS, (sum(rT^2) - tT * sum(diag(Koo))) / length(lab),
                         (sum(rS[lab] * rT) - tX * sum(diag(Koo))) / length(lab),
                         (sum(rS^2) - tS * sum(diag(GRM))) / n_all,
                         var_T = mean(rT^2), var_S = mean(rS^2),
                         dbar_T = mean(diag(Koo)), dbar_S = mean(diag(GRM)))
}

rows <- list(); chk <- numeric(0)
for (k in seq_len(n_rep)) {
  x  <- rnorm(n_all); e <- draw_noise(coh, vc)
  Y  <- sqrt(pve_x * V_T) * x + e[, 1]
  Sv <- sqrt(pve_x * V_S) * x + e[, 2]
  Ga <- gene_drop_cohort(coh, st$K, maf)                  # alpha SNPs, blockwise to bound memory
  Sv <- Sv + alpha * rowSums(Ga)
  G0 <- gene_drop_cohort(coh, n_null, maf)                # null SNPs
  lab <- label_set(coh, st$miss, "pairs")
  mydf <- make_mydf(coh, Y, Sv, x, lab)
  ind <- one_per_family(coh)
  syn <- add_fastscan_cache(SynSurrG_typeIerror_step1(mydf))
  obs <- ObsG_typeIerror_step1(mydf)
  ssr <- add_fastscan_cache(SynSurr_typeIerror_step1(mydf, ind))
  score <- function(G, type) {
    out <- NULL
    for (b0 in seq(1L, ncol(G), by = 500L)) {
      Gb <- G[, b0:min(ncol(G), b0 + 499L), drop = FALSE]
      out <- rbind(out, cbind(
        SynSurrG = score_test_SynSurrG_fast(Gb, syn)$negative_log10_pval_SynSurrG,
        ObsG     = score_test_ObsG_multiply(Gb, obs)$negative_log10_pval_ObsG,
        SynSurr  = score_test_SynSurrG_fast(Gb, ssr, subset = ind, label = "SynSurr")$negative_log10_pval_SynSurr))
    }
    data.frame(rep = k, type = type, snp = seq_len(ncol(G)), out)
  }
  ra <- score(Ga, "alpha"); r0 <- score(G0, "null")
  ## per-SNP refit of the variance components for the first alpha SNPs
  W0 <- cbind(1, x)
  ra$SynSurrG_refit <- NA_real_
  for (j in seq_len(min(n_refit, st$K))) {
    v  <- mom_vc(mydf$Y_obs, Sv, cbind(W0, Ga[, j]))
    pj <- add_fastscan_cache(step1_fixed(mydf, vc_fixed = v))
    ra$SynSurrG_refit[j] <- score_test_SynSurrG_fast(Ga[, j, drop = FALSE], pj)$negative_log10_pval_SynSurrG
    if (k == 1 && j <= 3) {                               # invariance check with fixed components
      m2 <- mydf; m2$S <- Sv - alpha * Ga[, j]
      p2 <- add_fastscan_cache(step1_fixed(m2, vc_fixed = as.list(syn$vc$used)))
      chk <- c(chk, abs(score_test_SynSurrG_fast(Ga[, j, drop = FALSE], p2)$T_score_SynSurrG -
                        score_test_SynSurrG_fast(Ga[, j, drop = FALSE], syn)$T_score_SynSurrG))
    }
  }
  rows[[k]] <- rbind(ra, transform(r0, SynSurrG_refit = NA_real_))
  rows[[k]]$tau_T2 <- syn$vc$used[["tau_T2"]]; rows[[k]]$tau_TS <- syn$vc$used[["tau_TS"]]
  rows[[k]]$tau_S2 <- syn$vc$used[["tau_S2"]]
  cat(sprintf("rep %d done (tau_S2 %.3f)\n", k, syn$vc$used[["tau_S2"]]))
}
out <- do.call(rbind, rows)
attr(out, "settings") <- c(as.list(st), list(alpha = alpha, rho = vc$rho, n_all = n_all, seed = seed,
                           array = array_id, invariance_max_abs_diff = if (length(chk)) max(chk) else NA))
saveRDS(out, out_file)
cat("invariance check max |diff| =", if (length(chk)) max(chk) else NA, "\nsaved", out_file, "\n")

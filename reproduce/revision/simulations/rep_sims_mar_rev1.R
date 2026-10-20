## Follow-up of the S5 missingness-at-random design of rep_sims_rev1.R
## (Reviewer 2, comment 1: may the labelling depend on the surrogates?).
##   S5W : labelling depends on a variable W observed for everyone, not in the model,
##         and correlated with Y beyond the synthetic phenotype Yhat (= S here):
##           d = residual of the regression of Y on Yhat (all individuals, true Y),
##           W = r_w * d / sd(d) + sqrt(1 - r_w^2) * N(0, 1),
##         so that corr(W, Y | Yhat) = r_w (0.3, 0.5); P(labelled) = plogis(a + b z),
##         z = standardised W, b = 1, a solved so the labelled fraction is 1 - miss
##         (label_set(type = "mar_s") of sim_dgp_rev1.R with W in place of S).
##         "mar_wy" (settings 7-10): W = Y + kappa * N(0, 1), a noisy measurement of the target
##         itself, kappa = sd(d) * sqrt(1 / r_w^2 - 1) so that again corr(W, Y | Yhat) = r_w.
##         "mar_w" carries almost no genetic signal beyond Yhat (W depends on G only through
##         the (1 - c) beta_G G part of d, c = cov(Y, Yhat) / var(Yhat) ~ 0.82), "mar_wy" all of it.
##         MCAR settings are rerun with new seeds as a Monte Carlo check of Table S5.
##         Seed 20261016 + 4e6 + 1e3 * setting + array.
##   S5V : the six S5 settings regenerated with the S5 seeds (20261016 + 3e6 + 1e3 * setting
##         + array), so the data are identical to results_rev1/S5/; Syn-PALM and SynSurr are
##         refitted with the TRUE variance components in place of the method-of-moments ones.
## Methods (both studies): SynSurrG (Syn-PALM, MoM VCs), ObsG, OracleG, SynSurr, Obs as in
## rep_sims_rev1.R, plus SynSurrG_trueVC and SynSurr_trueVC (true VCs plugged in).
##   Rscript rep_sims_mar_rev1.R <S5W|S5V> <setting> <array> [n_rep=25]
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R")
                   source("simulation_code/rev1/sim_dgp_rev1.R") })
options(synsurrg.verbose = FALSE)
args     <- commandArgs(trailingOnly = TRUE)
study    <- args[1]
setting  <- as.integer(args[2])
array_id <- as.integer(args[3])
n_rep    <- if (length(args) >= 4) as.integer(args[4]) else 25L
maf <- 0.25; n_obs0 <- 3e4
stopifnot(study %in% c("S5W", "S5V"))

## ---- settings -----------------------------------------------------------------
S <- switch(study,
  S5W = rbind(
    data.frame(label = rep(c("random", "mar_w", "mar_w"), 2), r_w = rep(c(0, 0.3, 0.5), 2),
               miss = rep(c(0.5, 0.9), each = 3), rho_tot = 0.5, stringsAsFactors = FALSE),
    ## settings 7-10 (added after 1-6 had run): W a noisy measurement of Y itself
    data.frame(label = "mar_wy", r_w = rep(c(0.3, 0.5), 2), miss = rep(c(0.5, 0.9), each = 2),
               rho_tot = 0.5, stringsAsFactors = FALSE)),
  ## identical to S5 of rep_sims_rev1.R (same order, hence same setting index and seed)
  S5V = expand.grid(label = c("random", "mar_x", "mar_s"), miss = c(0.5, 0.9), rho_tot = 0.5,
                    stringsAsFactors = FALSE))
st <- S[setting, , drop = FALSE]
miss <- st$miss; lab_type <- st$label
vc <- vc_from_total(st$rho_tot)
pve_x <- 0.10; her <- 5e-4; n_snp <- 20L; n_null <- 200L
n_all <- round(n_obs0 / (1 - miss))
coh <- make_cohort(n_all, prop_related = 1, r = 0.5)

OUT <- sprintf("simulation_code/results_rev1/S5/%s/", c(S5W = "mar_w", S5V = "true_vc")[[study]])
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("%ssetting_%02d_array_%03d.rds", OUT, setting, array_id)
if (file.exists(out_file)) quit(save = "no")
seed <- 20261016L + 1000000L * c(S5V = 3L, S5W = 4L)[[study]] + 1000L * setting + array_id
set.seed(seed)
cat(study, "setting", setting, ":", paste(names(st), unlist(st), collapse = ", "),
    sprintf("| rho_g %.3f rho_e %.3f rho %.3f | n %d | seed %d\n", vc$rho_g, vc$rho_e, vc$rho, n_all, seed))

## ---- true variance components plugged into the step-1 fits -----------------------
## The step-1 functions of funcs_v2.R are used unchanged except that the estimated
## components are overwritten right after they are computed (inserted statement).
insert_after <- function(f, pattern, expr) {
  b <- as.list(body(f)); i <- which(vapply(b, function(s) grepl(pattern, paste(deparse(s), collapse = " ")), TRUE))
  stopifnot(length(i) == 1L)
  body(f) <- as.call(c(b[1:i], list(expr), b[-(1:i)])); f
}
true_used <- c(tau_T2 = vc$Gm[1, 1], tau_TS = vc$Gm[1, 2], tau_S2 = vc$Gm[2, 2],
               sigma_T2 = vc$Em[1, 1], sigma_TS = vc$Em[1, 2], sigma_S2 = vc$Em[2, 2])
SynSurrG_step1_trueVC <- insert_after(SynSurrG_typeIerror_step1, "^vc <- \\.vc\\$vc$", bquote({
  hat_tau_T2 <- .(true_used[["tau_T2"]]); hat_tau_TS <- .(true_used[["tau_TS"]])
  hat_tau_S2 <- .(true_used[["tau_S2"]]); hat_sigma_T2 <- .(true_used[["sigma_T2"]])
  hat_sigma_TS <- .(true_used[["sigma_TS"]]); hat_sigma_S2 <- .(true_used[["sigma_S2"]])
  vc$used <- .(true_used) }))
## SynSurr (one individual per family, unrelated): marginal (co)variances Gm + Em
marg <- vc$Gm + vc$Em
SynSurr_step1_trueVC <- insert_after(SynSurr_typeIerror_step1, "^hat_sigma_TS <- cov", bquote({
  hat_sigma_T2 <- .(marg[1, 1]); hat_sigma_S2 <- .(marg[2, 2]); hat_sigma_TS <- .(marg[1, 2]) }))

V_T <- vc$Gm[1, 1] + vc$Em[1, 1]
beta_G <- sqrt(her * V_T / (2 * maf * (1 - maf)))
beta_X <- sqrt(pve_x * V_T)
rows <- list()
for (k in seq_len(n_rep)) {
  ## same draws, in the same order, as rep_sims_rev1.R (S5V reproduces results_rev1/S5/)
  G  <- gene_drop_cohort(coh, n_snp + n_null, maf)
  x  <- rnorm(n_all)
  e  <- draw_noise(coh, vc)
  b  <- c(rep(beta_G, n_snp), rep(0, n_null))
  Y  <- as.numeric(G %*% b) + beta_X * x + e[, 1]
  Sv <- as.numeric(G %*% b) + beta_X * x + e[, 2]
  if (lab_type %in% c("mar_w", "mar_wy")) {
    d <- residuals(lm(Y ~ Sv))
    W <- if (lab_type == "mar_w") st$r_w * d / sd(d) + sqrt(1 - st$r_w^2) * rnorm(n_all)
         else Y + sd(d) * sqrt(1 / st$r_w^2 - 1) * rnorm(n_all)
    lab <- label_set(coh, miss, "mar_s", s = W)          # logit in standardised W, b = 1
    pc_W   <- cor(residuals(lm(W ~ Sv)), d)               # corr(W, Y | Yhat)
    pc_Wx  <- cor(residuals(lm(W ~ Sv + x)), residuals(lm(Y ~ Sv + x)))
  } else {
    lab <- label_set(coh, miss, lab_type, x = x, s = Sv)
    pc_W <- pc_Wx <- NA_real_
  }
  mydf <- make_mydf(coh, Y, Sv, x, lab)
  ind <- one_per_family(coh)
  syn  <- add_fastscan_cache(SynSurrG_typeIerror_step1(mydf))
  synT <- add_fastscan_cache(SynSurrG_step1_trueVC(mydf))
  ssr  <- add_fastscan_cache(SynSurr_typeIerror_step1(mydf, ind))
  ssrT <- add_fastscan_cache(SynSurr_step1_trueVC(mydf, ind))
  r <- list(
    SynSurrG        = score_test_SynSurrG_fast(G, syn),
    SynSurrG_trueVC = score_test_SynSurrG_fast(G, synT, label = "SynSurrG_trueVC"),
    ObsG            = score_test_ObsG_multiply(G, ObsG_typeIerror_step1(mydf)),
    OracleG         = score_test_OracleG_multiply(G, OracleG_typeIerror_step1(mydf)),
    SynSurr         = score_test_SynSurrG_fast(G, ssr, subset = ind, label = "SynSurr"),
    SynSurr_trueVC  = score_test_SynSurrG_fast(G, ssrT, subset = ind, label = "SynSurr_trueVC"),
    Obs             = score_test_Obs_multiply(G, Obs_typeIerror_step1(mydf, ind), ind))
  mom <- syn$vc$used
  for (m in names(r)) rows[[length(rows) + 1]] <- data.frame(
    rep = k, method = m, snp = seq_along(b), beta = b,
    hat_beta = r[[m]][[paste0("hat_beta_", m)]], var_hat_beta = r[[m]][[paste0("var_hat_beta_", m)]],
    nlp = r[[m]][[paste0("negative_log10_pval_", m)]], n_lab = length(lab),
    tau_T2 = mom[["tau_T2"]], tau_TS = mom[["tau_TS"]], tau_S2 = mom[["tau_S2"]],
    sigma_T2 = mom[["sigma_T2"]], sigma_TS = mom[["sigma_TS"]], sigma_S2 = mom[["sigma_S2"]],
    pcor_W_Y_given_Yhat = pc_W, pcor_W_Y_given_Yhat_x = pc_Wx)
}
out <- do.call(rbind, rows)
attr(out, "settings") <- c(as.list(st), list(study = study, rho = vc$rho, rho_g = vc$rho_g, rho_e = vc$rho_e,
                           Gm = vc$Gm, Em = vc$Em, true_vc = true_used, her = her, beta_G = beta_G,
                           pve_x = pve_x, n_all = n_all, n_pair = coh$n_pair, seed = seed, array = array_id))
saveRDS(out, out_file)
cat("saved", out_file, "\n")

## Revision (NC R1, reviewer: are the method-of-moments variance components adequate?).
## REML refit of the Syn-PALM variance components for ONE protein, then the Syn-PALM and
## Obs-GWAS scan of the protein's SNP subset with those components plugged in.
## Data, transforms, GRM, covariates and per-SNP computation are exactly those of
## code/synpalm_main_code/step2_protein_fast.R (and of its subset copy
## code/post_cf3/sens/step2_sens.R): step-1 cache <p>_synthetic_phenotype_cf_noPRS_noPC.rds,
## crossfit_folds_K5.rds, sparse GRM (x2), INT(Y) on the labelled, rank-INT(Yhat) on all,
## v2 array genotypes with NA -> 0 and fix_constant_columns, add_fastscan_cache +
## score_test_SynSurrG_fast, score_test_ObsG_multiply. Only the six variance components
## change: instead of the MoM estimates of SynSurrG_ablation_estimate / ObsG_ablation_estimate
## they are REML estimates with the same sparse GRM, from the same three fits as reml() of
## the VC simulations (simulation_code/rev1/variance_components_rev1.R):
##   tau_T2, sigma_T2  univariate REML of INT(Y) on the labelled
##   tau_S2, sigma_S2  univariate REML of rank-INT(Yhat) on all individuals
##   G, E (2x2)        bivariate REML of (INT(Y), rank-INT(Yhat)) on the labelled (G, E PSD)
## each with all covariates of the scan. REML is computed exactly through the block-wise
## eigen-decomposition of the sparse GRM (code/post_cf3/vcrefit/reml_spectral.R); GMMAT::glmmkin,
## the implementation of the simulations, is run as a check of the univariate fits for a few
## proteins (VCREFIT_GMMAT=1). The combined 2x2 blocks go through the same
## constrain_vc_bivariate() as production (a no-op when they are PSD; the flags are kept).
## Plugged in: tau_T2, tau_TS, sigma_T2, sigma_TS from the bivariate REML on the labelled
## (G, E constrained PSD), tau_S2, sigma_S2 from the univariate REML of Yhat on all. (The
## simulations' combination, T variances from the univariate fit, is stored as
## reml_raw_simspec; the unconstrained bivariate REML as reml_biv_free.)
## Obs-GWAS uses the univariate REML (tau_T2, sigma_T2) via constrain_vc_univariate().
##   Rscript code/post_cf3/vcrefit/vcrefit.R <protein_index>
## env VCREFIT_SOURCE  reml (default) | mom: mom plugs in the production MoM estimates
##                     through the same code path (validation: must reproduce the main scan)
##     VCREFIT_SNP_FILE  rsids to scan (default results/vcrefit_cf_noPRS_noPC/snps/snps_<p>.txt)
## Outputs (results/vcrefit_cf_noPRS_noPC/):
##   vc/<p>_reml.rds                 REML estimates, convergence, timings, MoM for comparison
##   scan/<p>_<source>_scan_subset.rds  Syn-PALM / Obs-GWAS -log10 P, beta, var(beta)
## Finished steps are skipped on restart (preempted serial_requeue tasks just rerun).
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(Matrix); library(GMMAT) })
setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
S1_TAG  <- "_cf_noPRS_noPC"
SOURCE  <- Sys.getenv("VCREFIT_SOURCE", "reml")
stopifnot(SOURCE %in% c("reml", "mom"))
t_start <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")
protein_index <- as.integer(commandArgs(TRUE)[1])
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
OUT <- Sys.getenv("VCREFIT_OUT", "results/vcrefit_cf_noPRS_noPC/")   # test runs: another directory
dir.create(paste0(OUT, "vc"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(OUT, "scan"), recursive = TRUE, showWarnings = FALSE)
protein_names <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
this_protein  <- protein_names[protein_index]
vc_file   <- sprintf("%svc/%s_reml.rds", OUT, this_protein)
scan_file <- sprintf("%sscan/%s_%s_scan_subset.rds", OUT, this_protein, SOURCE)
snp_file  <- Sys.getenv("VCREFIT_SNP_FILE", sprintf("%ssnps/snps_%s.txt", OUT, this_protein))
if (file.exists(scan_file)) { tick("already done: ", scan_file); quit(save = "no") }

## ---------------------------------------------------------------- data (as step2_protein_fast.R)
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
df_allUKB_eur <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE,
                       data.table = FALSE, select = c("eid", "UKBcentre"))
df_all_variables_eur_orignal <- fread("data/extracted_pheno_data528_eur.tsv",
                                      sep = "\t", header = TRUE, data.table = FALSE)
FD  <- readRDS("data/crossfit_folds_K5.rds")
idx <- match(FD$analysis_eid, df_all_variables_eur_orignal$eid)
df_all_variables_eur <- df_all_variables_eur_orignal[idx, ]
rm(df_all_variables_eur_orignal); invisible(gc())
s1_file <- paste0("results/allukb_protein_pQTL_european_discovery/", this_protein, "/",
                  this_protein, "_synthetic_phenotype", S1_TAG, ".rds")
S1 <- readRDS(s1_file)
stopifnot(identical(as.character(S1$eid), as.character(FD$analysis_eid)),
          !isTRUE(S1$permuted), identical(S1$use_prs, FALSE), identical(S1$use_pcs, FALSE))
cat(sprintf("protein %s (index %d): rho_oof = %.4f, n_obs = %d, VC source = %s\n",
            this_protein, protein_index, S1$rho_oof, S1$n_obs, SOURCE))
y_all <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE, data.table = FALSE,
               select = c("eid", this_protein))
y_all <- y_all[match(FD$analysis_eid, y_all$eid), this_protein]
GRM <- sgrm[FD$index_all, FD$index_all]
rm(sgrm); invisible(gc())
cov_id  <- c("p21022", paste0("p22009_a", 1:20), "p31", "p22000")
test_df <- df_all_variables_eur[, c("eid", cov_id)]
colnames(test_df) <- c("eid", "age", paste0("geneticPC", 1:20), "sex", "batch")
test_df <- test_df %>% left_join(df_allUKB_eur, by = "eid")
test_df$age2     <- test_df$age^2
test_df$age_sex  <- test_df$age * (as.numeric(test_df$sex) - 1)
test_df$age2_sex <- test_df$age2 * (as.numeric(test_df$sex) - 1)
test_df$protein     <- y_all
test_df$protein_hat <- S1$protein_hat
stopifnot(all(as.character(test_df$eid) == as.character(FD$analysis_eid)))
rm(df_all_variables_eur, df_allUKB_eur); invisible(gc())
test_df <- INT(test_df, "protein")
n <- nrow(test_df)
r <- rank(test_df %>% pull(protein_hat))
test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
test_df_dmm     <- model.matrix(~ sex + batch + UKBcentre, data = test_df)
continuous_vars <- test_df[, setdiff(colnames(test_df), c("sex", "batch", "UKBcentre"))]
test_df_transformed <- cbind(continuous_vars, test_df_dmm)
test_df_transformed$`(Intercept)` <- NULL
colnames(test_df_transformed) <- gsub("-", "_", colnames(test_df_transformed))
obs_protein_index <- which(!is.na(test_df_transformed$protein))
X_all_names <- setdiff(colnames(test_df_transformed),
                       c("eid", "protein", "protein_hat", "int", "yhat_int"))
mydf <- list(X_all = test_df_transformed[, X_all_names],
             S     = test_df_transformed$yhat_int,
             Y_obs = test_df_transformed$int,
             GRM   = GRM)
tick("data ready (n = ", n, ", n_obs = ", length(obs_protein_index), ", ", length(X_all_names), " covariates)")

## MoM estimates exactly as production (only the six moments; cheap)
mom_raw <- function() {
  X  <- cbind(1, as.matrix(mydf$X_all))
  o  <- obs_protein_index
  rT <- as.numeric(lm.fit(X[o, ], mydf$Y_obs[o])$residuals)
  rS <- as.numeric(lm.fit(X, mydf$S)$residuals)
  Goo <- GRM[o, o]; Go <- Goo; diag(Go) <- 0; Ga <- GRM; diag(Ga) <- 0
  tT <- as.numeric(t(rT) %*% Go %*% rT) / sum(Go^2)
  tS <- as.numeric(t(rS) %*% Ga %*% rS) / sum(Ga^2)
  tX <- as.numeric(t(rS[o]) %*% Go %*% rT) / sum(Go^2)
  list(raw = c(tau_T2 = tT, tau_TS = tX, tau_S2 = tS,
               sigma_T2 = (sum(rT^2) - tT * sum(diag(Goo))) / length(o),
               sigma_TS = (sum(rS[o] * rT) - tX * sum(diag(Goo))) / length(o),
               sigma_S2 = (sum(rS^2) - tS * sum(diag(GRM))) / n),
       var_T = mean(rT^2), var_S = mean(rS^2), dbar_T = mean(diag(Goo)), dbar_S = mean(diag(GRM)))
}

## ---------------------------------------------------------------- REML (exact, spectral; code/post_cf3/vcrefit/reml_spectral.R)
## Same three fits as reml() of simulation_code/rev1/variance_components_rev1.R:
## univariate Y (labelled), univariate Yhat (all), bivariate (Y, Yhat) (labelled) -> tau_TS, sigma_TS.
## GMMAT::glmmkin is the implementation used there; its univariate fits are run as a
## check when VCREFIT_GMMAT=1 (Yhat on all individuals takes ~20 min); its bivariate fit
## returned exactly zero cross-covariances here (boundary refit) and is not used.
if (!file.exists(vc_file)) {
  source("code/post_cf3/vcrefit/reml_spectral.R")
  memb <- integer(n)
  for (b in seq_along(FD$blocks)) memb[FD$blocks[[b]]] <- b
  stopifnot(all(memb > 0))
  Xc <- cbind(1, as.matrix(mydf$X_all))
  ## aliased covariates (e.g. a batch absent among the labelled) are dropped, as lm() does
  nalias <- function(rows) { q <- qr(Xc[rows, ]); sort(q$pivot[seq_len(q$rank)]) }
  o  <- obs_protein_index
  ko <- nalias(o); ka <- nalias(seq_len(n))
  t0 <- Sys.time()
  So <- spec_decomp(GRM, o, memb)
  Sa <- spec_decomp(GRM, seq_len(n), memb)
  tick("eigen-decomposition of the GRM blocks: labelled min/max eigenvalue ", signif(min(So$lam), 3), "/", signif(max(So$lam), 3),
       ", all ", signif(min(Sa$lam), 3), "/", signif(max(Sa$lam), 3))
  rT <- reml_uni(So, mydf$Y_obs[o], Xc[o, ko])
  rSo <- reml_uni(So, mydf$S[o], Xc[o, ko])
  rS <- reml_uni(Sa, mydf$S, Xc[, ka])
  tick("REML Y (labelled): tau = ", signif(rT$tau, 4), ", sigma = ", signif(rT$sigma, 4),
       " | Yhat (all): tau = ", signif(rS$tau, 4), ", sigma = ", signif(rS$sigma, 4),
       " | Yhat (labelled): tau = ", signif(rSo$tau, 4), ", sigma = ", signif(rSo$sigma, 4))
  ## start of the bivariate fit: univariate fits, correlations from the MoM residual covariances
  M <- mom_raw()
  cg <- max(min(M$raw[["tau_TS"]], 0.9 * sqrt(rT$tau * rSo$tau)), -0.9 * sqrt(rT$tau * rSo$tau))
  ce <- max(min(M$raw[["sigma_TS"]], 0.9 * sqrt(rT$sigma * rSo$sigma)), -0.9 * sqrt(rT$sigma * rSo$sigma))
  st <- list(G = matrix(c(rT$tau + 1e-3, cg, cg, rSo$tau + 1e-3), 2), E = matrix(c(rT$sigma, ce, ce, rSo$sigma), 2))
  r2 <- reml_biv(So, NULL, NULL, NULL, start = st, rot_y = rT$rotated$ys, rot_s = rSo$rotated$ys, rot_X = rT$rotated$Xs)
  tick("REML bivariate (labelled): G = ", paste(signif(r2$G[c(1, 2, 4)], 4), collapse = " "),
       "; E = ", paste(signif(r2$E[c(1, 2, 4)], 4), collapse = " "), "; convergence ", r2$convergence,
       "; nll of the two starts ", paste(sprintf("%.3f", r2$nll_starts), collapse = " "))
  ## unconstrained bivariate REML (entries free): does REML also leave the PSD cone?
  r2u <- reml_biv(So, NULL, NULL, NULL, start = list(G = r2$G, E = r2$E), rot_y = rT$rotated$ys,
                  rot_s = rSo$rotated$ys, rot_X = rT$rotated$Xs, psd = FALSE)
  tick("REML bivariate unconstrained: G = ", paste(signif(r2u$G[c(1, 2, 4)], 4), collapse = " "),
       "; E = ", paste(signif(r2u$E[c(1, 2, 4)], 4), collapse = " "), "; nll ", sprintf("%.3f", r2u$nll),
       " (PSD fit ", sprintf("%.3f", r2$nll), ")")
  ## VCs plugged into Syn-PALM: (T, TS) blocks of the bivariate REML (PSD), S from all individuals
  raw <- c(tau_T2 = r2$G[1, 1], tau_TS = r2$G[1, 2], tau_S2 = rS$tau,
           sigma_T2 = r2$E[1, 1], sigma_TS = r2$E[1, 2], sigma_S2 = rS$sigma)
  ## as reml() of the VC simulations (T variances from the univariate fit), kept for reference
  raw_sim <- c(tau_T2 = rT$tau, tau_TS = r2$G[1, 2], tau_S2 = rS$tau,
               sigma_T2 = rT$sigma, sigma_TS = r2$E[1, 2], sigma_S2 = rS$sigma)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  gm <- NULL
  if (identical(Sys.getenv("VCREFIT_GMMAT", "0"), "1")) {
    ids <- paste0("i", seq_len(n)); K <- GRM; dimnames(K) <- list(ids, ids)
    K <- as(forceSymmetric(K), "CsparseMatrix")
    d <- data.frame(id = ids, y = mydf$Y_obs, s = mydf$S, Xc[, -1], check.names = TRUE)
    xn <- colnames(d)[-(1:3)]
    gT <- glmmkin(reformulate(xn[ko[-1] - 1], response = "y"), data = d[o, ], kins = K[o, o], id = "id", family = gaussian(), verbose = FALSE)
    gS <- glmmkin(reformulate(xn[ka[-1] - 1], response = "s"), data = d, kins = K, id = "id", family = gaussian(), verbose = FALSE)
    gm <- c(tau_T2 = gT$theta[[2]], sigma_T2 = gT$theta[[1]], tau_S2 = gS$theta[[2]], sigma_S2 = gS$theta[[1]])
    tick("GMMAT check: ", paste(sprintf("%s=%.5f", names(gm), gm), collapse = " "))
    rm(K, d); invisible(gc())
  }
  saveRDS(list(protein = this_protein, index = protein_index, rho_oof = S1$rho_oof, n = n, n_obs = length(o),
               reml_raw = raw, reml_raw_simspec = raw_sim, reml_uni_T = c(tau = rT$tau, sigma = rT$sigma),
               reml_biv = list(G = r2$G, E = r2$E, convergence = r2$convergence, nll = r2$nll, nll_starts = r2$nll_starts),
               reml_biv_free = list(G = r2u$G, E = r2u$E, convergence = r2u$convergence, nll = r2u$nll),
               reml_uni_S_labelled = c(tau = rSo$tau, sigma = rSo$sigma),
               boundary = c(T = rT$boundary, S = rS$boundary),
               converged = c(T = TRUE, S = TRUE, TS = r2$convergence == 0), secs = secs, gmmat_check = gm,
               mom_raw = M$raw, var_T = M$var_T, var_S = M$var_S, dbar_T = M$dbar_T, dbar_S = M$dbar_S,
               lam_range = c(labelled = range(So$lam), all = range(Sa$lam)),
               method = paste("exact REML via block-wise eigen-decomposition of the sparse GRM (x2):",
                              "univariate T (labelled), S (all); bivariate (T, S) on the labelled with G, E PSD -> TS")),
          vc_file)
  rm(So, Sa, rT, rS, rSo); invisible(gc())
  tick("REML saved: ", vc_file)
}
V <- readRDS(vc_file)
cat("MoM raw : ", paste(sprintf("%s=%+.4f", names(V$mom_raw), V$mom_raw), collapse = " "), "\n")
cat("REML raw: ", paste(sprintf("%s=%+.4f", names(V$reml_raw), V$reml_raw), collapse = " "), "\n")

## ---------------------------------------------------------------- null models with given VCs
## SynSurrG_ablation_estimate() of code/funcs_v2.R from "constrain the variance components"
## on, with the six raw estimates supplied instead of computed.
SynSurrG_fixedvc <- function(mydf, raw, var_T, var_S, dbar_T, dbar_S) {
  GRM <- mydf$GRM; n <- nrow(GRM)
  obs_protein_index <- which(!is.na(mydf$Y_obs))
  GRM_obs_obs  <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index, 1:n]
  n_obs <- length(obs_protein_index)
  .vc <- constrain_vc_bivariate(raw[["tau_T2"]], raw[["tau_TS"]], raw[["tau_S2"]],
                                raw[["sigma_T2"]], raw[["sigma_TS"]], raw[["sigma_S2"]],
                                var_T = var_T, var_S = var_S, label = paste0("fixedvc_", SOURCE),
                                dbar_T = dbar_T, dbar_S = dbar_S)
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  I_matrix <- sparseMatrix(i = seq_len(n_obs), j = obs_protein_index,
                           x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n))
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  inv_Sigma22 <- matrix_inv_block(wait_matrix = Sigma22)
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  Y <- mydf$Y_obs[obs_protein_index]
  hatY <- mydf$S
  X_all <- as.matrix(data.frame(intercept = rep(1, n), mydf$X_all))
  X_obs <- X_all[obs_protein_index, , drop = FALSE]
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att", "A22", "Atb", "B1", "B_mat",
              "B2_1", "B2", "bt2", "Btt1", "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12", "vc")
  P <- mget(params, envir = environment())
  vcu <- vc$used
  add_fastscan_cache(P, Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(n))
}
## ObsG_ablation_estimate() of code/funcs_v2.R with (tau_T2, sigma_T2) supplied
ObsG_fixedvc <- function(mydf, tau2, sigma2, var_T, dbar_T) {
  GRM <- mydf$GRM
  obs_protein_index <- which(!is.na(mydf$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  n_obs <- length(obs_protein_index)
  .vc <- constrain_vc_univariate(tau2, sigma2, var_T = var_T, label = paste0("ObsG_fixedvc_", SOURCE), dbar = dbar_T)
  vc <- .vc$vc
  Sigma11 <- .vc$tau2 * GRM_obs_obs + .vc$sigma2 * Diagonal(n_obs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)
  Y <- mydf$Y_obs[obs_protein_index]
  X <- cbind(1, as.matrix(mydf$X_all[obs_protein_index, , drop = FALSE]))
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  mget(c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X", "vc"), envir = environment())
}
raw_use <- if (SOURCE == "reml") V$reml_raw else V$mom_raw
## Obs-GWAS: univariate MoM of Y is the same (tau_T2, sigma_T2) as in the bivariate MoM
obs_use <- if (SOURCE == "reml") c(tau_T2 = V$reml_uni_T[["tau"]], sigma_T2 = V$reml_uni_T[["sigma"]]) else raw_use[c("tau_T2", "sigma_T2")]
Psyn <- SynSurrG_fixedvc(mydf, raw_use, V$var_T, V$var_S, V$dbar_T, V$dbar_S)
tick("Syn-PALM null model (", SOURCE, "): used ", paste(sprintf("%s=%+.4f", names(Psyn$vc$used), Psyn$vc$used), collapse = " "),
     "; projected_G = ", Psyn$vc$projected_G, ", projected_E = ", Psyn$vc$projected_E)
Pobs <- ObsG_fixedvc(mydf, obs_use[["tau_T2"]], obs_use[["sigma_T2"]], V$var_T, V$dbar_T)
tick("Obs-GWAS null model")
rm(GRM, mydf); invisible(gc())

## ---------------------------------------------------------------- scan of the SNP subset (as step2_sens.R)
G0 <- BEDMatrix::BEDMatrix(path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
rs_index <- intersect(colnames(G0), readLines(snp_file))
matched_indices <- match(test_df_transformed$eid, rownames(G0))
if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
fix_idx <- intersect(obs_protein_index, FD$independent_indices)
tick("scanning ", length(rs_index), " SNPs of ", snp_file, " on ", NCORES, " cores")
cn <- c("SynSurrG_negative_log10_pvalues", "SynSurrG_hat_beta", "SynSurrG_var_hat_beta",
        "ObsG_negative_log10_pvalues", "ObsG_hat_beta", "ObsG_var_hat_beta")
scan_block <- function(idx) {
  Gm <- as.matrix(G0[matched_indices, idx])
  Gm <- replace(Gm, is.na(Gm), 0)
  Gm <- fix_constant_columns(Gm, fix_idx)
  r1 <- score_test_SynSurrG_fast(Gm, Psyn, label = "SynSurrG")
  r2 <- score_test_ObsG_multiply(Gm, Pobs)
  out <- cbind(r1$negative_log10_pval_SynSurrG, r1$hat_beta_SynSurrG, r1$var_hat_beta_SynSurrG,
               r2$negative_log10_pval_ObsG, r2$hat_beta_ObsG, r2$var_hat_beta_ObsG)
  dimnames(out) <- list(colnames(Gm), cn)
  out
}
part_file <- sub("\\.rds$", "_partial.rds", scan_file)
blocks <- split(rs_index, ceiling(seq_along(rs_index) / 200))
done <- if (file.exists(part_file)) readRDS(part_file) else list()
todo <- setdiff(seq_along(blocks), as.integer(names(done)))
for (grp in split(todo, ceiling(seq_along(todo) / (4 * NCORES)))) {
  res <- parallel::mclapply(blocks[grp], function(b) tryCatch(scan_block(b), error = function(e) e),
                            mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) !is.matrix(x), logical(1))
  if (any(bad)) stop(sum(bad), " block(s) failed; first: ",
                     if (inherits(res[[which(bad)[1]]], "error")) conditionMessage(res[[which(bad)[1]]]) else "worker killed?")
  for (k in seq_along(grp)) done[[as.character(grp[k])]] <- res[[k]]
  saveRDS(done, part_file)
  tick(length(done), "/", length(blocks), " blocks")
}
results <- do.call(rbind, done[as.character(seq_along(blocks))])
stopifnot(identical(rownames(results), rs_index))
attr(results, "rho_oof") <- S1$rho_oof
attr(results, "n_obs")   <- S1$n_obs
attr(results, "vc_source") <- SOURCE
attr(results, "vc")      <- Psyn$vc
attr(results, "vc_obs")  <- Pobs$vc
attr(results, "code")    <- "code/post_cf3/vcrefit/vcrefit.R (2026-10-20)"
saveRDS(results, scan_file)
invisible(file.remove(part_file))
tick("Protein ", this_protein, " done: ", scan_file)

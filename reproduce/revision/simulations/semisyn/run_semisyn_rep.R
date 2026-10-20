## =====================================================================
## run_semisyn_rep.R -- ONE replicate of the semi-synthetic simulation
## (design fixed by prep_semisyn.R; see its header).
##
## Per replicate (seed 20261100 + rep, stored in the output):
##   Y*  = h(S0) + sqrt(V_G) g + sum_causal b_j G_j + e,     Var(Y*) ~ 1
##         h: fixed function of 4 real surrogates (BMI, SBP, HbA1c, LDL; rank-
##            normalised): linear + BMI x HbA1c product + SBP hinge, Var(h) = V_H
##         g ~ N(0, GRM) on the sparse GRM blocks, e ~ N(0, V_E)
##   S'  = S0 with planted effects: S_null SNPs on their column, causal_YS SNPs
##         and sqrt(share) * g (standardised) on Z = platelet count. Y* uses the
##         unmodified S0, so S_null SNPs move S' (and Yhat) but not Y*.
##   90% of Y* masked at random; the full step-1 pipeline is then run on S'
##   (and the ICD-10 table) inside the UKB-PPP relatedness-cluster folds:
##   screening of 100 continuous surrogates by |r| + 3 fixed predictors +
##   10 ICD-10 categories by Wilcoxon p, ranger (NTREE trees), out-of-fold Yhat;
##   then the step-2 model (INT of Y_obs and Yhat, covariates age, age^2, sex,
##   age x sex, age^2 x sex, 20 PCs, batch, centre) with Syn-PALM (SynSurrG) and
##   Obs-GWAS (ObsG), and a naive OLS GWAS of Yhat in everyone (contrast only),
##   on the fixed SNP set (design.rds).
## Output: simulation_code/results_rev1/semisyn/reps/rep_<rep>.rds
##   (per-SNP matrix + attributes: seed, rho_oof, VC, INT slope kappa, timings)
##
##   Rscript run_semisyn_rep.R <rep> [NTREE=300]
## =====================================================================
setwd("<PROJECT_DIR>")
suppressMessages({library(data.table); library(dplyr); library(Matrix); library(ranger); library(withr)})
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
options(synsurrg.verbose = FALSE)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

args  <- commandArgs(trailingOnly = TRUE)
rep_id <- as.integer(args[1])
NTREE <- if (length(args) >= 2) as.integer(args[2]) else 300L
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
OUT   <- "simulation_code/results_rev1/semisyn/"
RDIR  <- paste0(OUT, Sys.getenv("SEMISYN_SUBDIR", "reps"), "/")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("%srep_%03d.rds", RDIR, rep_id)
if (file.exists(out_file)) { cat("already done:", out_file, "\n"); quit(save = "no") }

## variance components of Y* (sum = 1) and the other fixed settings
V_H <- as.numeric(Sys.getenv("SEMISYN_VH", "0.23"))      # h(S0)
V_G <- 0.20                                              # polygenic (relatedness) term
V_C <- 0.06                                              # planted causal SNPs (fixed by design)
V_E <- 1 - V_H - V_G - V_C
Z_G_SHARE <- 0.30          # share of Var(S'_Z) contributed by g
MISS <- 0.9; K <- 5
SEED <- 20261100L + rep_id

design <- readRDS(paste0(OUT, "design.rds"))
DI     <- readRDS(paste0(OUT, "design_indiv.rds"))
snps   <- design$snps
FD     <- readRDS("data/crossfit_folds_K5.rds")
stopifnot(all(as.numeric(DI$eid) == as.numeric(FD$analysis_eid)))
n <- length(FD$analysis_eid); fold_vec <- FD$fold

### surrogate tables in the analysis order, rows re-assigned (prep step 1) ###
Sraw <- fread("data/extracted_pheno_data528_eur.tsv", sep = "\t", header = TRUE, data.table = FALSE)
Draw <- fread("data/disease_data_level1_eur.tsv", sep = "\t", header = TRUE, data.table = FALSE)
stopifnot(all(Draw$eid == Sraw$eid))
idx  <- match(FD$analysis_eid, Sraw$eid); stopifnot(!anyNA(idx))
Sraw <- Sraw[idx, ]; Draw <- Draw[idx, ]
## covariates for the GWAS model: REAL rows (they belong with the genotypes)
cov_id  <- c("p21022", paste0("p22009_a", 1:20), "p31", "p22000")
test_df <- Sraw[, c("eid", cov_id)]
colnames(test_df) <- c("eid", "age", paste0("geneticPC", 1:20), "sex", "batch")
cen <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE, data.table = FALSE, select = c("eid", "UKBcentre"))
test_df <- test_df %>% left_join(cen, by = "eid")
stopifnot(all(as.numeric(test_df$eid) == as.numeric(FD$analysis_eid)))
## predictors: re-assigned rows (S0), eid column kept in the analysis order
S0 <- Sraw[DI$perm, ]; S0$eid <- Sraw$eid
D0 <- Draw[DI$perm, ]; D0$eid <- Draw$eid
stopifnot(all(S0$p31 == Sraw$p31))
rm(Sraw, Draw, cen); invisible(gc())
tick("tables ready, n = ", n)

### GRM and g ~ N(0, GRM) blockwise ##########################################
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
GRM <- sgrm[FD$index_all, FD$index_all]
stopifnot(all(colnames(GRM) == unname(as.character(FD$analysis_eid))))
rm(sgrm); invisible(gc())
rgrm <- function(GRM, blocks) {          # one draw of N(0, GRM), negative eigenvalues set to 0
  z <- rnorm(nrow(GRM)); g <- numeric(nrow(GRM))
  sz <- lengths(blocks); bid <- rep.int(seq_along(blocks), sz)[order(unlist(blocks))]
  sm <- summary(as(GRM, "generalMatrix")); sm <- sm[bid[sm$i] == bid[sm$j], ]
  stopifnot(nrow(summary(as(GRM, "generalMatrix"))) == nrow(sm))   # GRM is block diagonal in `blocks`
  one <- unlist(blocks[sz == 1]); dg <- diag(GRM)
  g[one] <- sqrt(pmax(dg[one], 0)) * z[one]
  sm <- sm[sz[bid[sm$i]] > 1, ]; si <- sm$i; sj <- sm$j; sx <- sm$x
  for (ix in split(seq_along(si), bid[si])) {
    b <- blocks[[bid[si[ix[1]]]]]; b1 <- b[1] - 1L
    stopifnot(all(diff(b) == 1L))
    A <- matrix(0, length(b), length(b)); A[cbind(si[ix] - b1, sj[ix] - b1)] <- sx[ix]
    e <- eigen(A, symmetric = TRUE)
    g[b] <- e$vectors %*% (sqrt(pmax(e$values, 0)) * crossprod(e$vectors, z[b]))
  }
  g
}

### phenotypes for this replicate ##########################################
set.seed(SEED)
zt <- function(x) { r <- rank(x); qnorm((r - 0.375) / (length(x) - 0.75 + 1)) }
H <- design$H_COLS
zb <- zt(S0[[H["bmi"]]]); zs <- zt(S0[[H["sbp"]]]); za <- zt(S0[[H["hba1c"]]]); zl <- zt(S0[[H["ldl"]]])
h_raw <- 1.0 * zb + 0.6 * zs + 0.5 * za + 0.5 * zl + 0.5 * zb * za + 0.8 * pmax(zs - 0.5, 0)
h <- sqrt(V_H) * (h_raw - mean(h_raw)) / sd(h_raw)
g <- rgrm(GRM, FD$blocks); g <- (g - mean(g)) / sd(g)
tick("g drawn")
e <- rnorm(n, 0, sqrt(V_E))
Ystar <- h + sqrt(V_G) * g + DI$score[, "Y"] + e
## observed surrogates S': planted S-acting effects; g enters Z (platelets)
Sp <- S0
for (cc in names(design$H_COLS)) Sp[[design$H_COLS[cc]]] <- S0[[design$H_COLS[cc]]] + DI$score[, cc]
zc <- design$Z_COL[["plt"]]
Sp[[zc]] <- S0[[zc]] + DI$score[, "plt"] + design$sdS0[["plt"]] * sqrt(Z_G_SHARE / (1 - Z_G_SHARE)) * g
obs_index <- sort(sample.int(n, round((1 - MISS) * n)))
y_all <- rep(NA_real_, n); y_all[obs_index] <- Ystar[obs_index]
var_parts <- c(h = var(h), g = V_G * var(g), planted = var(DI$score[, "Y"]), e = var(e), Y = var(Ystar))
tick("Y* drawn; var parts: ", paste(sprintf("%s=%.3f", names(var_parts), var_parts), collapse = " "),
     "; labelled ", length(obs_index))

### step 1: cross-fitted RF exactly as step1_protein_crossfit.R ##############
cand <- setdiff(colnames(Sp)[8:ncol(Sp)], "prs")
cand <- cand[!grepl("^p22009_", cand)]
X_obs <- as.matrix(Sp[obs_index, cand])
Z_obs <- as.matrix(D0[obs_index, -1])
y_obs <- y_all[obs_index]
fixed_pred <- c("p20116_i0", "p21000_i0", "p31")
sanitize <- function(nm) gsub("[- ]", "_", nm)
screen_fold <- function(train_rows) {
  tp <- match(train_rows, obs_index); stopifnot(!anyNA(tp))
  r <- suppressWarnings(cor(X_obs[tp, , drop = FALSE], y_obs[tp], use = "pairwise.complete.obs"))
  r[is.na(r)] <- 0
  cov_pred <- c(colnames(X_obs)[order(abs(r), decreasing = TRUE)[1:100]], fixed_pred)
  yv <- y_obs[tp]
  sig <- vapply(seq_len(ncol(Z_obs)), function(i) {
    zi <- Z_obs[tp, i]; vals <- unique(zi[!is.na(zi)])
    if (length(vals) == 2) -log10(wilcox.test(x = yv[zi == vals[2]], y = yv[zi == vals[1]],
                                               alternative = "two.sided", exact = FALSE)$p.value) else 0
  }, numeric(1))
  list(cov_pred = cov_pred, cov_dis = colnames(Z_obs)[order(sig, decreasing = TRUE)[1:10]])
}
build_frame <- function(rows, cov_pred, cov_dis) {
  a <- Sp[rows, cov_pred, drop = FALSE]
  b <- D0[rows, cov_dis, drop = FALSE]
  b[] <- lapply(b, function(x) factor(x, levels = c(0, 1), labels = c("0", "1")))
  out <- cbind(a, b); names(out) <- sanitize(names(out)); out
}
yhat <- rep(NA_real_, n); feats <- vector("list", K)
for (k in seq_len(K)) {
  train_rows <- intersect(obs_index, which(fold_vec != k))
  pred_rows  <- which(fold_vec == k)
  sel <- screen_fold(train_rows); feats[[k]] <- sel
  tr <- build_frame(train_rows, sel$cov_pred, sel$cov_dis); tr$protein <- y_all[train_rows]
  keep <- vapply(tr[setdiff(names(tr), "protein")], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  if (any(!keep)) tr <- tr[, c(setdiff(names(tr), c(names(keep)[!keep], "protein")), "protein"), drop = FALSE]
  rf <- ranger(protein ~ ., data = tr, num.trees = NTREE, save.memory = FALSE, write.forest = TRUE,
               seed = SEED * 10L + k, num.threads = NCORES)
  pr <- build_frame(pred_rows, sel$cov_pred, sel$cov_dis)[, setdiff(names(tr), "protein"), drop = FALSE]
  yhat[pred_rows] <- predict(rf, data = pr, num.threads = NCORES)$predictions
  rm(tr, pr, rf); invisible(gc())
}
stopifnot(!anyNA(yhat))
rho_oof <- cor(y_obs, yhat[obs_index])
rho_all <- cor(Ystar, yhat)
fold_rho <- sapply(seq_len(K), function(k) { ii <- intersect(obs_index, which(fold_vec == k)); cor(y_all[ii], yhat[ii]) })
top_feats <- sort(table(unlist(lapply(feats, function(f) f$cov_pred[1:10]))), decreasing = TRUE)
tick("cross-fit done: rho_oof = ", round(rho_oof, 4), " (all = ", round(rho_all, 4), "), by fold ",
     paste(round(fold_rho, 3), collapse = " "), "; top features: ", paste(names(top_feats)[1:8], collapse = ","))

### INT slope: d INT(y)/dy averaged over the labelled Y* ######################
yo <- Ystar[obs_index]; dens <- density(yo, n = 4096)
fy <- approx(dens$x, dens$y, xout = yo)$y
ty <- qnorm((rank(yo) - 0.375) / (length(yo) - 0.75 + 1))
kappa <- mean(fy / dnorm(ty))

### step 2 model frame exactly as step2_protein_fast.R #######################
test_df$age2 <- test_df$age^2
test_df$age_sex  <- test_df$age * (as.numeric(test_df$sex) - 1)
test_df$age2_sex <- test_df$age2 * (as.numeric(test_df$sex) - 1)
test_df$protein <- y_all; test_df$protein_hat <- yhat
test_df <- INT(test_df, "protein")
r <- rank(test_df$protein_hat); test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
## batch (p22000) and centre enter as their numeric codes, exactly as in the
## production step2_protein_fast.R (kept identical on purpose, 2026-10-20)
test_df_dmm <- model.matrix(~ sex + batch + UKBcentre, data = test_df)
continuous_vars <- test_df[, setdiff(colnames(test_df), c("sex", "batch", "UKBcentre"))]
tdt <- cbind(continuous_vars, test_df_dmm); tdt$`(Intercept)` <- NULL
colnames(tdt) <- gsub("-", "_", colnames(tdt))
obs_protein_index <- which(!is.na(tdt$protein)); stopifnot(identical(obs_protein_index, obs_index))
X_all_names <- setdiff(colnames(tdt), c("eid", "protein", "protein_hat", "int", "yhat_int"))
mydf <- list(X_all = tdt[, X_all_names], S = tdt$yhat_int, Y_obs = tdt$int, GRM = GRM)
pSG <- SynSurrG_ablation_estimate(mydf); vcu <- pSG$vc$used
pSG <- add_fastscan_cache(pSG, Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(n))
pOG <- ObsG_ablation_estimate(mydf)
## naive OLS GWAS of INT(Yhat) in everyone (contrast: what a Yhat-only GWAS sees)
Xn <- cbind(1, as.matrix(tdt[, X_all_names])); qrX <- qr(Xn)
ry <- qr.resid(qrX, tdt$yhat_int); df_res <- n - qrX$rank
tick("null models fitted; VC used: ", paste(sprintf("%s=%.3f", names(vcu), unlist(vcu)), collapse = " "))
vc_save <- list(SynSurrG = pSG$vc, ObsG = pOG$vc)
rm(GRM, mydf, continuous_vars, test_df, test_df_dmm, S0, Sp, D0, X_obs, Z_obs); invisible(gc())

### scan of the fixed SNP set ##################################################
G0 <- BEDMatrix::BEDMatrix("genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
matched <- match(unname(as.character(FD$analysis_eid)), rownames(G0)); stopifnot(!anyNA(matched))
fix_idx <- intersect(obs_protein_index, FD$independent_indices)
cn <- c("SynPALM_mlog10p", "SynPALM_beta", "SynPALM_var", "ObsG_mlog10p", "ObsG_beta", "ObsG_var",
        "YhatOLS_mlog10p", "YhatOLS_beta", "YhatOLS_var")
scan_block <- function(ids) {
  Gm <- as.matrix(G0[matched, ids, drop = FALSE]); Gm <- replace(Gm, is.na(Gm), 0)
  Gm <- with_seed(SEED + 7L, fix_constant_columns(Gm, fix_idx))
  a <- score_test_SynSurrG_fast(Gm, pSG, label = "SynSurrG")
  b <- score_test_ObsG_multiply(Gm, pOG)
  Gr <- qr.resid(qrX, Gm); ss <- colSums(Gr^2); bo <- colSums(Gr * ry) / ss
  s2 <- (sum(ry^2) - bo^2 * ss) / (df_res - 1); vo <- s2 / ss
  po <- -(pchisq(bo^2 / vo, 1, lower.tail = FALSE, log.p = TRUE)) / log(10)
  out <- cbind(a$negative_log10_pval_SynSurrG, a$hat_beta_SynSurrG, a$var_hat_beta_SynSurrG,
               b$negative_log10_pval_ObsG, b$hat_beta_ObsG, b$var_hat_beta_ObsG, po, bo, vo)
  dimnames(out) <- list(ids, cn); out
}
blocks_snp <- split(snps$rsid, ceiling(seq_along(snps$rsid) / 200))
res <- parallel::mclapply(blocks_snp, function(i) tryCatch(scan_block(i), error = function(e) e),
                          mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(res, function(x) !is.matrix(x), logical(1))
if (any(bad)) stop(sum(bad), " block(s) failed; first: ", conditionMessage(res[[which(bad)[1]]]))
res <- do.call(rbind, res); stopifnot(identical(rownames(res), snps$rsid))
attr(res, "settings") <- list(rep = rep_id, seed = SEED, ntree = NTREE, V_H = V_H, V_G = V_G, V_C = V_C,
                              V_E = V_E, Z_G_SHARE = Z_G_SHARE, miss = MISS, n = n, n_obs = length(obs_index),
                              h = "1.0 zBMI + 0.6 zSBP + 0.5 zHbA1c + 0.5 zLDL + 0.5 zBMI*zHbA1c + 0.8 (zSBP-0.5)_+",
                              design_seeds = design$seeds)
attr(res, "rho") <- c(oof = rho_oof, all = rho_all, by_fold = fold_rho)
attr(res, "kappa") <- kappa
attr(res, "sd_Yobs") <- sd(yo)
attr(res, "var_parts") <- var_parts
attr(res, "vc") <- vc_save
attr(res, "top_features") <- top_feats
attr(res, "time_sec") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
attr(res, "code") <- "run_semisyn_rep.R (2026-10-20)"
saveRDS(res, out_file)
tick("saved ", out_file)

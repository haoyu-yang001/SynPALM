## REAL-DATA POP-GWAS AND LABELLED-ONLY ANALYSIS (revision 2026-10-20; reply to
## Reviewer 1, Comments 1 and 3, and Reviewer 3, minor 1). Copy of
## code/post_cf3/sens/step2_sens.R (itself a copy of code/synpalm_main_code/
## step2_protein_fast.R): same data, covariates, inverse-normal transformations,
## cross-fitted Yhat (step-1 cache _cf_noPRS_noPC), GRM and Syn-PALM / Obs-GWAS score
## tests, on the SNP set results/sens_cf_noPRS_noPC/snps_<protein>.txt. Modes
## (env REALPOP_MODE):
##   full             all 402,086 analysis individuals, as in the main analysis.
##                    Besides Syn-PALM and Obs-GWAS, the three GWAS that POP-GWAS
##                    combines are computed, each as a mixed model with the GRM
##                    restricted to its own sample (as in simulation_code/rev1/
##                    typeIerror_split_rev1.R) and as ordinary linear regression:
##                      Z1  Yhat on the labelled individuals
##                      Z2  Yhat on the unlabelled individuals
##                      Z3  Y on the labelled individuals (mixed = Obs-GWAS)
##                    POP-GWAS itself is computed from these in summarize_realpop.R.
##   full_jointPERM1  as full, after the joint permutation of (INT Y, INT Yhat) over
##                    all analysis individuals with seed 1, exactly as in
##                    code/synpalm_main_code/step2_protein_fast_jointperm.R (Fig S16):
##                    every SNP is null; the labelled set becomes a random subset.
##   labonly          the unlabelled individuals are dropped: the analysis sample
##                    is the labelled individuals only, GRM restricted to them, Yhat
##                    re-inverse-normal-transformed within them. Syn-PALM and
##                    Obs-GWAS only (Obs-GWAS must equal the main Obs-GWAS).
##   Rscript step2_realpop.R <protein_index>
## Output: results/realpop_cf_noPRS_noPC/<protein>_<mode>.rds (one row per SNP)
library(data.table)
library(dplyr)
library(Matrix)

setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
S1_TAG <- "_cf_noPRS_noPC"
MODE   <- Sys.getenv("REALPOP_MODE", "full")
stopifnot(MODE %in% c("full", "full_jointPERM1", "labonly"))
POP    <- MODE != "labonly"
t_start <- Sys.time()
tick <- function(...) cat(sprintf("[%6.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")

args          <- commandArgs(trailingOnly = TRUE)
protein_index <- as.integer(args[1])
NCORES        <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
OUTDIR        <- "results/realpop_cf_noPRS_noPC/"
dir.create(OUTDIR, showWarnings = FALSE)

sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])

df_allUKB_eur <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE,
                       data.table = FALSE, select = c("eid", "UKBcentre"))
df_all_variables_eur_orignal <- fread("data/extracted_pheno_data528_eur.tsv",
                                      sep = "\t", header = TRUE, data.table = FALSE)
protein_names <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
this_protein  <- protein_names[protein_index]
final_file <- paste0(OUTDIR, this_protein, "_", MODE, ".rds")
if (file.exists(final_file)) { tick("exists: ", final_file); quit(save = "no") }

FD <- readRDS("data/crossfit_folds_K5.rds")
idx <- match(FD$analysis_eid, df_all_variables_eur_orignal$eid)
df_all_variables_eur <- df_all_variables_eur_orignal[idx, ]
rm(df_all_variables_eur_orignal); gc()

s1_file <- paste0("results/allukb_protein_pQTL_european_discovery/", this_protein, "/",
                  this_protein, "_synthetic_phenotype", S1_TAG, ".rds")
S1 <- readRDS(s1_file)
stopifnot(identical(as.character(S1$eid), as.character(FD$analysis_eid)), isFALSE(S1$permuted),
          identical(S1$use_prs, FALSE), identical(S1$use_pcs, FALSE))
cat(sprintf("mode %s; step-1 cache %s: rho_oof = %.4f, n_obs = %d\n", MODE, basename(s1_file), S1$rho_oof, S1$n_obs))

y_all <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE, data.table = FALSE,
               select = c("eid", this_protein))
y_all <- y_all[match(FD$analysis_eid, y_all$eid), this_protein]

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
rm(df_all_variables_eur, df_allUKB_eur); gc()

test_df <- INT(test_df, "protein")
n <- nrow(test_df)
r <- rank(test_df %>% pull(protein_hat))
test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))

test_df_dmm     <- model.matrix(~ sex + batch + UKBcentre, data = test_df)
continuous_vars <- test_df[, setdiff(colnames(test_df), c("sex", "batch", "UKBcentre"))]
test_df_transformed <- cbind(continuous_vars, test_df_dmm)
test_df_transformed$`(Intercept)` <- NULL
colnames(test_df_transformed) <- gsub("-", "_", colnames(test_df_transformed))

if (MODE == "full_jointPERM1") {                 # as step2_protein_fast_jointperm.R
  set.seed(1)
  perm <- sample(nrow(test_df_transformed))
  test_df_transformed$int      <- test_df_transformed$int[perm]
  test_df_transformed$yhat_int <- test_df_transformed$yhat_int[perm]
  cat(sprintf("*** JOINT PERMUTATION of (Y, Yhat), seed 1: cor(Y, Yhat) among labeled = %.4f ***\n",
              cor(test_df_transformed$int, test_df_transformed$yhat_int, use = "complete.obs")))
}

GRM  <- sgrm[FD$index_all, FD$index_all]
rm(sgrm); invisible(gc())
keep <- seq_len(n)
independent_indices <- FD$independent_indices
if (MODE == "labonly") {
  keep <- which(!is.na(test_df_transformed$int))
  test_df_transformed <- test_df_transformed[keep, ]
  GRM  <- GRM[keep, keep]
  ## Yhat re-transformed within the labelled sample (Y's INT is over the labelled
  ## individuals already, so Y is unchanged)
  nk <- length(keep); rk <- rank(test_df_transformed$protein_hat)
  test_df_transformed$yhat_int <- qnorm((rk - 0.375) / (nk - 2 * 0.375 + 1))
  independent_indices <- which(keep %in% FD$independent_indices)
}
obs_protein_index <- which(!is.na(test_df_transformed$int))
unl_index         <- which(is.na(test_df_transformed$int))
X_all_names <- setdiff(colnames(test_df_transformed),
                       c("eid", "protein", "protein_hat", "int", "yhat_int"))
X_all <- test_df_transformed[, X_all_names]
## covariate columns constant within a sample are dropped from the models of that
## sample (main / full Syn-PALM and Obs-GWAS: none dropped, so they equal the main scan)
nonconst <- function(ix) X_all_names[vapply(X_all[ix, ], function(v) length(unique(v)) > 1, logical(1))]
if (MODE == "labonly") {
  kx <- nonconst(seq_len(nrow(X_all)))
  cat("labonly: dropped constant covariates:", paste(setdiff(X_all_names, kx), collapse = ", "), "\n")
  X_all <- X_all[, kx]
}
mydf <- list(X_all = X_all, S = test_df_transformed$yhat_int, Y_obs = test_df_transformed$int, GRM = GRM)
cat(sprintf("n = %d, labelled = %d, unlabelled = %d, covariates = %d\n",
            nrow(X_all), length(obs_protein_index), length(unl_index), ncol(X_all)))

tick("data ready")
SynSurrG_step1_pars <- SynSurrG_ablation_estimate(mydf)
vcu <- SynSurrG_step1_pars$vc$used
SynSurrG_step1_pars <- add_fastscan_cache(SynSurrG_step1_pars,
                         Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(nrow(GRM)))
tick("Syn-PALM null model")
ObsG_step1_pars <- ObsG_ablation_estimate(mydf)
tick("Obs-GWAS null model")

if (POP) {
  ## mixed-model GWAS of Yhat in the labelled (Z1) and the unlabelled (Z2) sample,
  ## GRM restricted to that sample: ObsG_ablation_estimate with Y_obs = Yhat there
  S <- test_df_transformed$yhat_int
  mk <- function(ix, cols) list(X_all = X_all[, cols], S = S, GRM = GRM,
                                Y_obs = replace(rep(NA_real_, length(S)), ix, S[ix]))
  cols_unl <- nonconst(unl_index)
  cat("unlabelled sample: dropped constant covariates:", paste(setdiff(X_all_names, cols_unl), collapse = ", "), "\n")
  Z1_pars <- ObsG_ablation_estimate(mk(obs_protein_index, X_all_names))
  Z2_pars <- ObsG_ablation_estimate(mk(unl_index, cols_unl))
  stopifnot(identical(Z1_pars$obs_protein_index, obs_protein_index),
            identical(Z2_pars$obs_protein_index, unl_index))
  tick("POP-GWAS mixed null models")
  ## ordinary least squares pieces (residual, sigma2, X, (X'X)^-1) for one sample
  ols <- function(ix, y, cols) {
    X <- cbind(1, as.matrix(X_all[ix, cols])); y <- y[ix]
    XtX_inv <- solve(crossprod(X)); res <- as.numeric(y - X %*% (XtX_inv %*% crossprod(X, y)))
    list(ix = ix, X = X, XtX_inv = XtX_inv, res = res, sigma2 = sum(res^2) / (length(y) - ncol(X)))
  }
  L_y1 <- ols(obs_protein_index, test_df_transformed$int, X_all_names)
  L_s1 <- ols(obs_protein_index, S, X_all_names)
  L_s2 <- ols(unl_index, S, cols_unl)
  r_pheno <- cor(L_y1$res, L_s1$res)
  cat(sprintf("POP-GWAS: corr(Y, Yhat | covariates) among labelled = %.4f\n", r_pheno))
  ols_test <- function(G, L) {
    g <- G[L$ix, , drop = FALSE]
    A <- crossprod(L$X, g)
    gmg <- colSums(g * g) - colSums(A * (L$XtX_inv %*% A))
    cbind(as.numeric(crossprod(g, L$res)) / gmg, L$sigma2 / gmg)
  }
}
rm(GRM, mydf); invisible(gc())

G0 <- BEDMatrix::BEDMatrix(path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
snp_file  <- sprintf("results/sens_cf_noPRS_noPC/snps_%s.txt", this_protein)
rs_index  <- intersect(colnames(G0), readLines(snp_file))
tick("scanning ", length(rs_index), " SNPs of ", snp_file)
chunk_size <- 200L
matched_indices <- match(test_df_transformed$eid, rownames(G0))
if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
fix_idx <- intersect(obs_protein_index, independent_indices)

scan_block <- function(idx) {
  Gm <- as.matrix(G0[matched_indices, idx])
  Gm <- replace(Gm, is.na(Gm), 0)
  Gm <- fix_constant_columns(Gm, fix_idx)
  r1 <- score_test_SynSurrG_fast(Gm, SynSurrG_step1_pars, label = "SynSurrG")
  r2 <- score_test_ObsG_multiply(Gm, ObsG_step1_pars)
  out <- cbind(SynSurrG_negative_log10_pvalues = r1$negative_log10_pval_SynSurrG,
               SynSurrG_hat_beta = r1$hat_beta_SynSurrG, SynSurrG_var_hat_beta = r1$var_hat_beta_SynSurrG,
               ObsG_negative_log10_pvalues = r2$negative_log10_pval_ObsG,
               ObsG_hat_beta = r2$hat_beta_ObsG, ObsG_var_hat_beta = r2$var_hat_beta_ObsG)
  if (POP) {
    z1 <- score_test_ObsG_multiply(Gm, Z1_pars); z2 <- score_test_ObsG_multiply(Gm, Z2_pars)
    out <- cbind(out,
                 Z1mm_beta = z1$hat_beta_ObsG, Z1mm_var = z1$var_hat_beta_ObsG,
                 Z2mm_beta = z2$hat_beta_ObsG, Z2mm_var = z2$var_hat_beta_ObsG,
                 `colnames<-`(ols_test(Gm, L_y1), c("Z3lm_beta", "Z3lm_var")),
                 `colnames<-`(ols_test(Gm, L_s1), c("Z1lm_beta", "Z1lm_var")),
                 `colnames<-`(ols_test(Gm, L_s2), c("Z2lm_beta", "Z2lm_var")),
                 eaf = colMeans(Gm) / 2)
  }
  rownames(out) <- colnames(Gm)
  out
}

blocks <- split(rs_index, ceiling(seq_along(rs_index) / chunk_size))
res <- parallel::mclapply(blocks, function(idx) tryCatch(scan_block(idx), error = function(e) e),
                          mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(res, function(x) !is.matrix(x), logical(1))
if (any(bad)) stop(sum(bad), " block(s) failed; first: ", conditionMessage(res[[which(bad)[1]]]))
results <- do.call(rbind, res)
stopifnot(identical(rownames(results), rs_index))
attr(results, "mode")    <- MODE
attr(results, "rho_oof") <- S1$rho_oof
attr(results, "n_total") <- nrow(test_df_transformed)
attr(results, "n_lab")   <- length(obs_protein_index)
attr(results, "n_unlab") <- length(unl_index)
attr(results, "vc")      <- SynSurrG_step1_pars$vc
attr(results, "vc_obs")  <- ObsG_step1_pars$vc
if (POP) {
  attr(results, "r_pheno") <- r_pheno
  attr(results, "vc_Z1") <- Z1_pars$vc; attr(results, "vc_Z2") <- Z2_pars$vc
}
attr(results, "code") <- "code/post_cf3/realpop/step2_realpop.R (2026-10-20)"
saveRDS(results, final_file)
tick("Protein ", this_protein, " ", MODE, " done: ", final_file)

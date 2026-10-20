## SENSITIVITY COPY (revision 2026-10-19; reply to Reviewer 1, Comments 4 and 8) of
## code/synpalm_main_code/step1_protein_crossfit.R: identical, except that
##   SENS_FOLD_FILE  fold assignment to use (default data/crossfit_folds_K5.rds)
##   SENS_NPRED      number of screened continuous predictors (default 100)
##   SENS_TAG        suffix appended to the run tag of every output (required)
## =====================================================================
## STEP 1 of 2 -- cross-fitted synthetic phenotype for ONE protein.
##
## Why this is a separate job: the fold loop now trains 5 random forests on
## ~80% of the labeled samples each, roughly 40x the work of the old single
## model on 10%. The old script repeated all of that inside every chunk_array
## job. Doing it once per protein and caching the result removes that waste.
##
##   Rscript step1_protein_crossfit.R <protein_index>
##
## Run make_folds_once.R first.
## =====================================================================

library(data.table)
library(dplyr)
library(Matrix)
library(ranger)

setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")

args          <- commandArgs(trailingOnly = TRUE)
protein_index <- as.integer(args[1])

## ---- switches -------------------------------------------------------  ## <<< CF
USE_PRS   <- FALSE   # reply letter, Comment 8: no polygenic score in the predictor set.
# The 2026-08-14 run had this switch but the PRS blocks below were
# hard-coded `if (TRUE)`, so the PRS entered every fold anyway
# (that version is kept in code/before/step1_protein_crossfit_PRSbug_2026-08-14.R).
USE_PCS   <- FALSE   # genetic PCs (p22009_a1-a20) are not RF predictors either;
# they stay in the GWAS covariates (step 2).
PERMUTE_Y <- identical(Sys.getenv("SYNPALM_PERMUTE_Y", "FALSE"), "TRUE")
#            TRUE = calibration check: permute the observed protein among
# labeled individuals BEFORE any model is fitted, then run the
# whole pipeline. Out-of-fold rho must come out ~ 0.
PERM_SEED <- as.integer(Sys.getenv("SYNPALM_PERM_SEED", "1"))
## every output of this run carries RUN_TAG, so files from the 2026-08 run
## (PRS in the forest) can never be picked up by mistake
RUN_TAG   <- paste0("_cf_noPRS_noPC", Sys.getenv("SENS_TAG"))
stopifnot(nzchar(Sys.getenv("SENS_TAG")))
SENS_NPRED <- as.integer(Sys.getenv("SENS_NPRED", "100"))
cat("sensitivity run", RUN_TAG, " n predictors", SENS_NPRED, " folds", Sys.getenv("SENS_FOLD_FILE", "data/crossfit_folds_K5.rds"), "\n")

## ranger otherwise grabs every core it can see, which oversubscribes a shared
## node when --cpus-per-task is small. Take the allocation SLURM gave us.     ## <<< CF
NTHREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(NTHREADS) || NTHREADS < 1) NTHREADS <- 1
cat("ranger threads:", NTHREADS, "\n")

cat("protein index:", protein_index, "\n")

## ---- cached fold assignment -----------------------------------------  ## <<< CF
## Step 1 needs no GRM: the fold assignment and the analysis-cohort ordering
## were fixed once by make_folds_once.R and are the same for every protein.
FD       <- readRDS(Sys.getenv("SENS_FOLD_FILE", "data/crossfit_folds_K5.rds"))
K        <- FD$K
fold_vec <- FD$fold                            # length n_all, values 1..K
n_all    <- length(fold_vec)

## ---- phenotype data --------------------------------------------------
## Only the single protein column is read, not the full 2,920-column table.  ## <<< CF
hdr           <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))
protein_names <- hdr[49:2968]
this_protein  <- protein_names[protein_index]

df_allUKB_eur_protein <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE,
                               data.table = FALSE, select = c("eid", this_protein))
df_all_variables_eur_orignal <- fread("data/extracted_pheno_data528_eur.tsv",
                                      sep = "\t", header = TRUE, data.table = FALSE)
disease_data_level1_eur_orignal <- fread("data/disease_data_level1_eur.tsv",
                                         sep = "\t", header = TRUE, data.table = FALSE)
stopifnot(all(disease_data_level1_eur_orignal$eid == df_all_variables_eur_orignal$eid))

idx <- match(FD$analysis_eid, df_all_variables_eur_orignal$eid)
stopifnot(!any(is.na(idx)))
df_all_variables_eur    <- df_all_variables_eur_orignal[idx, ]
disease_data_level1_eur <- disease_data_level1_eur_orignal[idx, ]
stopifnot(identical(as.character(df_allUKB_eur_protein$eid), as.character(FD$analysis_eid)),
          identical(as.character(df_all_variables_eur$eid),  as.character(FD$analysis_eid)),
          identical(as.vector(as.numeric(disease_data_level1_eur$eid)), as.vector(as.numeric(FD$analysis_eid))))

obs_index <- which(!is.na(df_allUKB_eur_protein[, this_protein]))   # labeled rows
cat("labeled individuals for", this_protein, ":", length(obs_index), "\n")

## ---- genotypes (only needed if the PRS is used) ----------------------  ## <<< CF
if (USE_PRS) {
  G0 <- BEDMatrix::BEDMatrix(
    path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
  geno_samp.id <- rownames(G0)
}

## ---- protein PRS (kept only to preserve column positions) ------------  ## <<< CF
G_PRS <- data.frame(prs = rep(0, nrow(df_allUKB_eur_protein)))
rownames(G_PRS) <- df_allUKB_eur_protein$eid
if (USE_PRS) {
  path_PRS <- "data/proteinPRSmodel"
  dirs <- tolower(list.dirs(path_PRS, full.names = FALSE, recursive = FALSE))
  if (this_protein %in% dirs) {
    matched_indices_prs <- match(df_allUKB_eur_protein$eid, geno_samp.id)
    target_path <- file.path(path_PRS, toupper(this_protein))
    txt_files   <- list.files(target_path, pattern = "\\.txt$", full.names = TRUE)
    combined_data <- do.call(rbind, lapply(txt_files, read.table,
                                           header = TRUE, sep = "\t", stringsAsFactors = FALSE))
    intercept       <- as.numeric(combined_data[1, 3])
    snp_weights_all <- combined_data[-1, ]
    common_snps     <- intersect(snp_weights_all[, 1], colnames(G0))
    matched_weights <- snp_weights_all[match(common_snps, snp_weights_all[, 1]), ]
    G0_selected     <- G0[matched_indices_prs, common_snps, drop = FALSE]
    G0_selected[is.na(G0_selected)] <- 0
    G_PRS <- G0_selected %*% as.numeric(matched_weights[, 3]) / 10 + intercept / 5
    colnames(G_PRS) <- "prs"
  }
}
rm(df_all_variables_eur_orignal, disease_data_level1_eur_orignal); gc()
stopifnot(all(df_all_variables_eur$eid == rownames(G_PRS)))
df_all_variables_eur <- data.frame(df_all_variables_eur, PRS = G_PRS)

## =====================================================================
## CROSS-FITTED PREDICTION                                               ## <<< CF
## The whole supervised pipeline -- covariate screening, ICD screening,
## random forest -- is refitted inside each fold. Screening uses the
## protein, so running it once on the pooled labeled data would let the
## fold-k model see fold-k phenotypes and would defeat cross-fitting.
## =====================================================================

y_all <- df_allUKB_eur_protein[, this_protein]

if (PERMUTE_Y) {                                                          ## <<< CF
  set.seed(PERM_SEED)
  y_all[obs_index] <- sample(y_all[obs_index])
  cat("*** PERMUTATION CHECK: protein permuted among labeled individuals ***\n")
}

## Screening only ever looks at labeled individuals, so the screening matrices
## are built over obs_index only. Materialising them over all ~400k rows would
## cost ~2 GB for nothing. Global row indices are mapped into these matrices
## with match(., obs_index).                                              ## <<< CF
## candidate continuous surrogates: columns 8.. of the predictor table, minus the
## PRS column added above and, unless USE_PCS, the 20 genetic PCs
cand <- setdiff(colnames(df_all_variables_eur)[8:ncol(df_all_variables_eur)], "prs")
if (!USE_PCS) cand <- cand[!grepl("^p22009_", cand)]
cat("candidate continuous surrogates:", length(cand), "\n")
X_obs <- as.matrix(df_all_variables_eur[obs_index, cand])
Z_obs <- as.matrix(disease_data_level1_eur[obs_index, -1])
y_obs <- y_all[obs_index]

## fixed predictors that are always kept (no screening involved)
fixed_pred <- c("p20116_i0", "p21000_i0", "p31")
if (USE_PRS) fixed_pred <- c(fixed_pred, "prs")

sanitize <- function(nm) gsub("[- ]", "_", nm)

screen_fold <- function(train_rows) {
  tp <- match(train_rows, obs_index)          # positions inside X_obs / Z_obs / y_obs
  stopifnot(!any(is.na(tp)))
  
  ## ---- continuous covariates: |Pearson r| with the protein -----------
  r <- suppressWarnings(cor(X_obs[tp, , drop = FALSE], y_obs[tp],
                            use = "pairwise.complete.obs"))
  r[is.na(r)] <- 0
  cov_pred <- c(colnames(X_obs)[order(abs(r), decreasing = TRUE)[seq_len(SENS_NPRED)]], fixed_pred)
  
  ## ---- ICD-10 level-1 categories: Wilcoxon -log10 p ------------------
  yv <- y_obs[tp]
  sig <- vapply(seq_len(ncol(Z_obs)), function(i) {
    zi   <- Z_obs[tp, i]
    vals <- unique(zi[!is.na(zi)])
    if (length(vals) == 2) {
      p <- wilcox.test(x = yv[zi == vals[2]], y = yv[zi == vals[1]],
                       alternative = "two.sided", paired = FALSE, exact = FALSE)$p.value
      -log10(p)
    } else 0
  }, numeric(1))
  
  ## NOTE ON THE ORIGINAL CODE: the old branch was
  ##   if (length(which(sig_icd10_log10 > 5) <= 5)) { top 5 } else { top 10 }
  ## The "<= 5" sits INSIDE which(), so length() just counts the significant
  ## categories and the if() treats any non-zero count as TRUE. The effect is
  ## the opposite of what it reads like: >=1 significant category -> top 5,
  ## none -> top 10. It also disagrees with the manuscript, which says 10.
  ## Written out as intended:
  ## Revision: always the 10 categories with the smallest p-values, as the
  ## Methods state (the 2026-08 run kept 5 when <= 5 categories had p < 1e-5).
  n_sig    <- sum(sig > 5)
  n_keep   <- 10L
  cov_dis  <- colnames(Z_obs)[order(sig, decreasing = TRUE)[seq_len(n_keep)]]
  
  list(cov_pred = cov_pred, cov_dis = cov_dis, n_sig = n_sig)
}

build_frame <- function(rows, cov_pred, cov_dis) {
  a <- df_all_variables_eur[rows, cov_pred, drop = FALSE]
  b <- disease_data_level1_eur[rows, cov_dis, drop = FALSE]
  b[] <- lapply(b, function(x) factor(x, levels = c(0, 1), labels = c("0", "1")))
  out <- cbind(a, b)
  names(out) <- sanitize(names(out))
  out
}

protein_hat <- rep(NA_real_, n_all)                                       ## <<< CF
fold_features <- vector("list", K)

for (k in seq_len(K)) {
  train_rows <- intersect(obs_index, which(fold_vec != k))   # labeled, outside fold k
  pred_rows  <- which(fold_vec == k)                         # EVERYONE in fold k
  cat(sprintf("fold %d: train on %d labeled, predict %d individuals\n",
              k, length(train_rows), length(pred_rows)))
  
  sel <- screen_fold(train_rows)
  fold_features[[k]] <- sel
  
  tr <- build_frame(train_rows, sel$cov_pred, sel$cov_dis)
  tr$protein <- y_all[train_rows]
  
  ## drop predictors that are constant or all-NA within this fold's training set,
  ## otherwise ranger fails or the factor levels do not match at prediction time
  keep <- vapply(tr[setdiff(names(tr), "protein")],
                 function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  drop_nm <- names(keep)[!keep]
  if (length(drop_nm)) {
    cat("  dropping non-varying predictors:", paste(drop_nm, collapse = ", "), "\n")
    tr <- tr[, c(setdiff(names(tr), c(drop_nm, "protein")), "protein"), drop = FALSE]
  }
  
  rf <- ranger(protein ~ ., data = tr, num.trees = 300, importance = "impurity",
               keep.inbag = TRUE, save.memory = FALSE, write.forest = TRUE,
               seed = 260611 + k, num.threads = NTHREADS)
  
  pr <- build_frame(pred_rows, sel$cov_pred, sel$cov_dis)
  pr <- pr[, setdiff(names(tr), "protein"), drop = FALSE]   # same columns, same order
  protein_hat[pred_rows] <- predict(rf, data = pr)$predictions
  
  rm(tr, pr, rf); gc()
}
stopifnot(!any(is.na(protein_hat)))

## ---- out-of-fold prediction accuracy ---------------------------------  ## <<< CF
cor_results <- cor(y_all[obs_index], protein_hat[obs_index])
fold_rho <- sapply(seq_len(K), function(k) {
  ii <- intersect(obs_index, which(fold_vec == k)); cor(y_all[ii], protein_hat[ii])
})
cat("out-of-fold rho =", round(cor_results, 4),
    " per fold:", paste(round(fold_rho, 4), collapse = " "), "\n")

outdir <- paste0("results/allukb_protein_pQTL_european_discovery/", this_protein, "/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tag    <- paste0(RUN_TAG, if (PERMUTE_Y) paste0("_PERMUTED", PERM_SEED) else "")
saveRDS(cor_results,        paste0(outdir, this_protein, "_correlation", tag, ".rds"))
saveRDS(fold_rho,           paste0(outdir, this_protein, "_correlation_by_fold", tag, ".rds"))
saveRDS(length(obs_index),  paste0(outdir, this_protein, "_number_of_observed_protein", tag, ".rds"))
saveRDS(fold_features,      paste0(outdir, this_protein, "_selected_features_by_fold", tag, ".rds"))

## ---- cache the synthetic phenotype for step 2 ------------------------  ## <<< CF
## Only a length-n numeric vector plus diagnostics: a few MB, not the
## step-1 parameter objects (V11^{-1} etc.), which are large.
saveRDS(list(eid           = df_allUKB_eur_protein$eid,
             protein_hat   = protein_hat,
             rho_oof       = cor_results,
             rho_by_fold   = fold_rho,
             n_obs         = length(obs_index),
             fold_features = fold_features,
             use_prs       = USE_PRS,
             use_pcs       = USE_PCS,
             permuted      = PERMUTE_Y,
             perm_seed     = if (PERMUTE_Y) PERM_SEED else NA,
             run_tag       = RUN_TAG),
        paste0(outdir, this_protein, "_synthetic_phenotype", tag, ".rds"))

if (PERMUTE_Y) {
  cat("PERMUTATION CHECK: out-of-fold rho =", round(cor_results, 4),
      "-- this must be close to 0 (|rho| < ~0.02 at n_obs ~ 33,000).\n")
} else {
  cat("step 1 done for", this_protein, " out-of-fold rho =", round(cor_results, 4), "\n")
}
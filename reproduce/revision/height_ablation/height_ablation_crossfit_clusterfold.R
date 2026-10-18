###############################################################################
## Syn-PALM ablation study -- UKB height, v2 genotypes (REVISION, 2026-10-13)
## Cross-fitted, ONE job per missing rate, folds by RELATEDNESS CLUSTER.
##
## Same cohort, masking seed, predictors (phenotypes only: no PRS, no genetic
## PCs) and RF settings as Synthetic_GWAS_ablation_v2geno_crossfit_prep.R; the
## changes are
##   * folds: the relatedness-cluster folds of data/crossfit_folds_K5.rds (the
##     UKB-PPP folds), so relatives are never split across folds; the 2026-09
##     run used sample(rep_len(1:K, n)), which split relative pairs.
##   * the genome-wide scan uses the vectorised score test
##     (code/fastscan/funcs_fastscan.R) inside this job, NSEG resumable segments.
##   * the oracle methods (OracleG, Oracle) do not depend on the missing rate and
##     are computed only in the missing_rate = 0 job.
## Output (same 18 column names as the 2026-09 run, Oracle columns NA for
## missing_rate > 0):
##   ablation_results/v2geno_main_results_crossfit_clusterfold/
##     height_missing_<r>_negative_log10_pvals_v2geno_crossfit_clusterfold.rds
##
## Usage:  Rscript height_ablation_crossfit_clusterfold.R <missing_rate> [NSEG=20]
###############################################################################
setwd("<PROJECT_DIR>/")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
t_start <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")
args <- commandArgs(trailingOnly = TRUE)
missing_rate <- as.numeric(args[1])
NSEG <- if (length(args) >= 2) as.integer(args[2]) else 20L
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
K <- 5; seed_miss <- 2025317; seed_rf <- 25313; seed_indep <- 25324; seed_fold_extra <- 20261013
out_dir <- "ablation_results/v2geno_main_results_crossfit_clusterfold"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
final_file <- file.path(out_dir, paste0("height_missing_", missing_rate,
                                        "_negative_log10_pvals_v2geno_crossfit_clusterfold.rds"))
if (file.exists(final_file)) { cat("already done:", final_file, "\n"); quit(save = "no") }
DO_ORACLE <- missing_rate == 0

### read the grm matrix ################
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
rm(sGRM)
sgrm <- sgrm*2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
#obtain the eid
eids_sGRM <- colnames(sgrm)
length(eids_sGRM)

G0 <- BEDMatrix::BEDMatrix(path = paste0("<PROJECT_DIR>/genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed",sep = ""), simple_names = T)

### read phenotype data ################
direction <- "<UKB_DIR>/Pharma_Proteomics_Project_2923_proteins_data/"
df_phenos <- fread(paste0(direction,"/phenotypes/extracted_pheno_data.tsv",sep=""),data.table = F )
colnames(df_phenos) <- gsub("^participant\\.", "", colnames(df_phenos))
#colnames(df_phenos)
df_phenos_colnames <- fread(paste0(direction,"/phenotypes/extracted_pheno_data_meaning.csv",sep=""),data.table = F )
colnames(df_phenos) <- df_phenos_colnames[-1,2]

###choose the european discovery sample and with the grm records
subeids0 <- df_phenos$eid[which(!is.na(df_phenos$`Genetic ethnic grouping`))]

subeids <- intersect(eids_sGRM, subeids0)
length(subeids) ## 408727
all(subeids %in% eids_sGRM)
all(subeids %in% df_phenos$eid)

### rank df_phenos as the eid order in grm subeids
df_phenos <- df_phenos %>% filter(eid %in% subeids)
df_phenos <- df_phenos[match(subeids, df_phenos$eid),]
all(df_phenos$eid == subeids)

df_phenos$Sex <- as.factor(df_phenos$Sex)
df_phenos$`Smoking status | Instance 0` <- factor(df_phenos$`Smoking status | Instance 0`,levels = c(0,1,2,-3))
df_phenos$`Time blood sample collected | Instance 0 | Array 0` <- NULL
## using the median to impute
for (i in 2:ncol(df_phenos)) {
  df_phenos[which(is.na(df_phenos[,i])),i] <- median(as.numeric(df_phenos[,i]),na.rm = T)
}

## continous first, then ordered
df_phenos0 <- data.frame(eid = df_phenos$eid,
                         age=df_phenos$`Age at recruitment`,
                         bmi=df_phenos$`Body mass index (BMI) | Instance 0`,
                         height = df_phenos$`Standing height | Instance 0`,
                         weight = df_phenos$`Weight | Instance 0`,
                         wbc = df_phenos$`Waist circumference | Instance 0`,
                         bfp = df_phenos$`Body fat percentage | Instance 0`,
                         chol = df_phenos$`Cholesterol | Instance 0`,
                         trigl = df_phenos$`Triglycerides | Instance 0`,
                         sbp = df_phenos$`Systolic blood pressure, automated reading | Instance 0 | Array 0`,
                         dbp = df_phenos$`Diastolic blood pressure, automated reading | Instance 0 | Array 0`,
                         HbA1c = df_phenos$`Glycated haemoglobin (HbA1c) | Instance 0`,
                         hdl = df_phenos$`HDL cholesterol | Instance 0`,
                         ldl = df_phenos$`LDL direct | Instance 0`,
                         wbcc = df_phenos$`White blood cell (leukocyte) count | Instance 0`,
                         rbc = df_phenos$`Red blood cell (erythrocyte) count | Instance 0`,
                         plet = df_phenos$`Platelet crit | Instance 0`,
                         lympt = df_phenos$`Lymphocyte percentage | Instance 0`,
                         monot = df_phenos$`Monocyte percentage | Instance 0`,
                         APOA = df_phenos$`Apolipoprotein A | Instance 0`,
                         APOB = df_phenos$`Apolipoprotein B | Instance 0`,
                         TDI = df_phenos$`Townsend deprivation index at recruitment`,
                         geneticPC1 = df_phenos$`Genetic principal components | Array 1`,
                         geneticPC2 = df_phenos$`Genetic principal components | Array 2`,
                         geneticPC3 = df_phenos$`Genetic principal components | Array 3`,
                         geneticPC4 = df_phenos$`Genetic principal components | Array 4`,
                         geneticPC5 = df_phenos$`Genetic principal components | Array 5`,
                         geneticPC6 = df_phenos$`Genetic principal components | Array 6`,
                         geneticPC7 = df_phenos$`Genetic principal components | Array 7`,
                         geneticPC8 = df_phenos$`Genetic principal components | Array 8`,
                         geneticPC9 = df_phenos$`Genetic principal components | Array 9`,
                         geneticPC10 = df_phenos$`Genetic principal components | Array 10`,
                         geneticPC11 = df_phenos$`Genetic principal components | Array 11`,
                         geneticPC12 = df_phenos$`Genetic principal components | Array 12`,
                         geneticPC13 = df_phenos$`Genetic principal components | Array 13`,
                         geneticPC14 = df_phenos$`Genetic principal components | Array 14`,
                         geneticPC15 = df_phenos$`Genetic principal components | Array 15`,
                         geneticPC16 = df_phenos$`Genetic principal components | Array 16`,
                         geneticPC17 = df_phenos$`Genetic principal components | Array 17`,
                         geneticPC18 = df_phenos$`Genetic principal components | Array 18`,
                         geneticPC19 = df_phenos$`Genetic principal components | Array 19`,
                         geneticPC20 = df_phenos$`Genetic principal components | Array 20`,
                         sex = df_phenos$Sex,
                         smoking = df_phenos$`Smoking status | Instance 0`,
                         alcohol = as.factor(df_phenos$`Alcohol intake frequency. | Instance 0`),
                         fast_time_collected = as.factor(df_phenos$`Fasting time | Instance 0`),
                         batch = as.factor(df_phenos$`Genotype measurement batch`),
                         UKBcentre = as.factor(df_phenos$`UK Biobank assessment centre | Instance 0`))

df_allUKB_eur0 <- df_phenos0[match(subeids, df_phenos0$eid),]
geno_eids <- rownames(G0)

df_allUKB_eur <- df_allUKB_eur0[df_allUKB_eur0$eid %in% geno_eids, ] # 401832     48


analysis_df <- df_allUKB_eur
rownames(analysis_df) <- NULL
n <- nrow(analysis_df)
cat("analysis set size n =", n, "\n")

### 1. mask the outcome FIRST (same seed as before) ##########################
draw <- with_seed(seed_miss, sample(seq_len(n), size = round(missing_rate * n), replace = FALSE))
analysis_df$height_obs <- analysis_df$height
if (length(draw) > 0) analysis_df$height_obs[draw] <- NA
n_obs_total <- sum(!is.na(analysis_df$height_obs))
cat("observed Y:", n_obs_total, "/", n, "\n")

### 2. relatedness-cluster folds ##############################################
index_all <- match(as.character(analysis_df$eid), eids_sGRM)
stopifnot(!any(is.na(index_all)))
GRM <- sgrm[index_all, index_all]
rm(sgrm); invisible(gc())
hh <- find_blocks_vectorized(GRM)
FD <- readRDS("data/crossfit_folds_K5.rds")
f0 <- FD$fold[match(as.character(analysis_df$eid), as.character(FD$analysis_eid))]
cat("individuals with a UKB-PPP fold:", sum(!is.na(f0)), "/", n, "\n")
fold <- integer(n)
set.seed(seed_fold_extra)
n_conflict <- 0L
for (b in hh) {
  fb <- unique(f0[b][!is.na(f0[b])])
  if (length(fb) > 1) n_conflict <- n_conflict + 1L
  fold[b] <- if (length(fb)) fb[1] else sample.int(K, 1)
}
stopifnot(all(fold %in% seq_len(K)))
cat("blocks:", length(hh), " blocks with >1 UKB-PPP fold:", n_conflict, "\n")
print(table(fold))
analysis_df$fold <- fold
stopifnot(all(vapply(hh, function(b) length(unique(fold[b])) == 1, logical(1))))

### 3. out-of-fold random forest predictions of height (unchanged) ###########
cov_pred <- c("age","sex","weight","wbc","bfp","smoking","chol","sbp","dbp","hdl","ldl")
analysis_df$height_hat <- NA_real_
for (k in seq_len(K)) {
  idx_pred  <- which(analysis_df$fold == k)
  idx_train <- which(analysis_df$fold != k & !is.na(analysis_df$height_obs))
  train_rf <- analysis_df[idx_train, cov_pred, drop = FALSE]
  train_rf$height <- analysis_df$height_obs[idx_train]
  rf_model <- ranger(height ~ ., data = train_rf, num.trees = 300,
                     save.memory = FALSE, write.forest = TRUE, seed = seed_rf + k,
                     num.threads = n_threads)
  analysis_df$height_hat[idx_pred] <-
    predict(rf_model, data = analysis_df[idx_pred, cov_pred, drop = FALSE],
            num.threads = n_threads)$predictions
  rm(rf_model, train_rf); invisible(gc())
}
stopifnot(!any(is.na(analysis_df$height_hat)))
obs_i <- !is.na(analysis_df$height_obs)
rho <- c(all = cor(analysis_df$height, analysis_df$height_hat),
         obs = if (sum(obs_i) > 2) cor(analysis_df$height[obs_i], analysis_df$height_hat[obs_i]) else NA,
         by_fold = sapply(seq_len(K), function(k) { i <- analysis_df$fold == k
                          cor(analysis_df$height[i], analysis_df$height_hat[i]) }))
tick("out-of-fold rho (all): ", round(rho[["all"]], 4))
saveRDS(rho, file.path(out_dir, paste0("height_miss_", missing_rate, "_correlation_clusterfold.rds")))

### 4. model frame (unchanged) ###############################################
cov_id <- c("age",paste0("geneticPC",1:20),"sex","batch","UKBcentre")
model_df <- analysis_df[, c("eid", cov_id)]
model_df$height     <- analysis_df$height
model_df$height_hat <- analysis_df$height_hat
model_df$height_obs <- analysis_df$height_obs
model_df <- INT(model_df, "height_obs")
r <- rank(model_df$height_hat); model_df$yhat_int    <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
r <- rank(model_df$height);     model_df$yoracle_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
model_df_dmm <- model.matrix(~ sex + batch + UKBcentre, data = model_df)
stopifnot(nrow(model_df_dmm) == n)
continuous_vars <- model_df[, setdiff(colnames(model_df), c("sex", "batch", "UKBcentre"))]
mt <- cbind(continuous_vars, model_df_dmm); mt$`(Intercept)` <- NULL
colnames(mt) <- gsub("-", "_", colnames(mt))
independent_indices <- with_seed(seed_indep, sapply(hh, function(x) x[sample.int(length(x), 1)]))
obs_protein_index <- which(!is.na(mt$height_obs))
X_all_names <- setdiff(colnames(mt), c("eid", "height", "height_hat", "height_obs", "int", "yhat_int", "yoracle_int"))
mydf <- list(X_all = mt[, X_all_names], S = mt$yhat_int, Y = mt$yoracle_int, Y_obs = mt$int, GRM = GRM)

### 5. null models ##############################################################
pSG <- SynSurrG_ablation_estimate(mydf); vcu <- pSG$vc$used
pSG <- add_fastscan_cache(pSG, Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(n))
pOG <- ObsG_ablation_estimate(mydf)
pS  <- add_fastscan_cache(SynSurr_ablation_estimate(mydf, independent_indices))
pO  <- Obs_ablation_estimate(mydf, independent_indices)
if (DO_ORACLE) {
  pRG <- OracleG_ablation_estimate(mydf)
  pR  <- Oracle_ablation_estimate(mydf, independent_indices)
}
tick("null models fitted")
vc_save <- list(SynSurrG = pSG$vc, ObsG = pOG$vc)
rm(GRM, mydf, model_df, continuous_vars, analysis_df, df_allUKB_eur); invisible(gc())

### 6. genome-wide scan ######################################################
matched_indices <- match(as.character(mt$eid), rownames(G0))
stopifnot(!any(is.na(matched_indices)))
fix_idx <- intersect(obs_protein_index, independent_indices)
cn <- c(t(outer(c("OracleG", "SynSurrG", "ObsG", "Oracle", "SynSurr", "Obs"),
                c("negative_log10_pvalues", "hat_beta", "var_hat_beta"), paste, sep = "_")))
rs_index <- colnames(G0)
seg_list <- split(seq_along(rs_index), cut(seq_along(rs_index), breaks = NSEG, labels = FALSE))
seg_file <- function(s) file.path(out_dir, paste0("height_missing_", missing_rate, "_seg", s, "of", NSEG, ".rds"))
NAc <- function(m) cbind(rep(NA_real_, m), rep(NA_real_, m), rep(NA_real_, m))
scan_block <- function(idx) {
  Gm <- as.matrix(G0[matched_indices, idx]); Gm <- replace(Gm, is.na(Gm), 0)
  Gm <- fix_constant_columns(Gm, fix_idx)
  a <- score_test_SynSurrG_fast(Gm, pSG, label = "SynSurrG")
  b <- score_test_ObsG_multiply(Gm, pOG)
  c <- score_test_SynSurrG_fast(Gm, pS, subset = independent_indices, label = "SynSurr")
  d <- score_test_Obs_multiply(Gm, pO, independent_indices)
  if (DO_ORACLE) {
    e <- score_test_OracleG_multiply(Gm, pRG); f <- score_test_Oracle_multiply(Gm, pR, independent_indices)
    eo <- cbind(e$negative_log10_pval_OracleG, e$hat_beta_OracleG, e$var_hat_beta_OracleG)
    fo <- cbind(f$negative_log10_pval_Oracle,  f$hat_beta_Oracle,  f$var_hat_beta_Oracle)
  } else { eo <- NAc(ncol(Gm)); fo <- NAc(ncol(Gm)) }
  out <- cbind(eo,
               a$negative_log10_pval_SynSurrG, a$hat_beta_SynSurrG, a$var_hat_beta_SynSurrG,
               b$negative_log10_pval_ObsG,     b$hat_beta_ObsG,     b$var_hat_beta_ObsG,
               fo,
               c$negative_log10_pval_SynSurr,  c$hat_beta_SynSurr,  c$var_hat_beta_SynSurr,
               d$negative_log10_pval_Obs,      d$hat_beta_Obs,      d$var_hat_beta_Obs)
  dimnames(out) <- list(colnames(Gm), cn); out
}
for (s in seq_len(NSEG)) {
  if (file.exists(seg_file(s))) { tick("segment ", s, " already done"); next }
  blocks <- split(seg_list[[s]], ceiling(seq_along(seg_list[[s]]) / 200))
  res <- parallel::mclapply(blocks, function(i) tryCatch(scan_block(i), error = function(e) e),
                            mc.cores = n_threads, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) !is.matrix(x), logical(1))
  if (any(bad)) stop("segment ", s, ": ", sum(bad), " block(s) failed; first: ", conditionMessage(res[[which(bad)[1]]]))
  saveRDS(do.call(rbind, res), seg_file(s))
  tick("segment ", s, "/", NSEG, " saved")
}
results <- do.call(rbind, lapply(seq_len(NSEG), function(s) readRDS(seg_file(s))))
stopifnot(identical(rownames(results), rs_index))
attr(results, "rho") <- rho; attr(results, "n") <- n; attr(results, "n_obs") <- length(obs_protein_index)
attr(results, "vc") <- vc_save; attr(results, "folds") <- "relatedness clusters (crossfit_folds_K5.rds)"
attr(results, "code") <- "height_ablation_crossfit_clusterfold.R (2026-10-13)"
saveRDS(results, final_file)
invisible(file.remove(vapply(seq_len(NSEG), seg_file, "")))
tick("done: ", final_file)

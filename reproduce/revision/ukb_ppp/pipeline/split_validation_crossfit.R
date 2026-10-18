## =====================================================================
## SPLIT-SAMPLE REPLICATION (revision; no PRS, no PCs in the forest).
##
## The analysis cohort is split 1:1 into a discovery and a validation half,
## by whole relatedness clusters (FD$blocks), so no pair of relatives is split
## across the halves. Inside each half the full step-1 pipeline is rerun on
## that half's data only: per-fold screening + random forest, 5-fold
## cross-fit with the global relatedness-cluster folds restricted to the half.
## Each half is then scanned genome-wide with Syn-PALM (SynSurrG) and the
## observed-only LMM (ObsG), with the same null models and fast score test
## as step2_protein_fast.R. The halves share no individual, so no protein
## value of one half reaches the synthetic phenotype of the other.
##
##   Rscript split_validation_crossfit.R <protein_index> [NSEG=10]
##
## Output: results/split_validation_cf_noPRS_noPC/<p>_split_scan.rds
##   591,558 x 12 matrix: {SynSurrG,ObsG} x {discovery,validation} x
##   {negative_log10_pvalues, hat_beta, var_hat_beta}
## =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(Matrix); library(ranger)
})
setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
t_start <- Sys.time()
tick <- function(...) cat(sprintf("[%6.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")

args          <- commandArgs(trailingOnly = TRUE)
protein_index <- as.integer(args[1])
NSEG          <- if (length(args) >= 2) as.integer(args[2]) else 10L
NCORES        <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
SEED_SPLIT    <- 250818L
RUN_TAG       <- "_split_cf_noPRS_noPC"
OUTDIR        <- "results/split_validation_cf_noPRS_noPC/"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

FD       <- readRDS("data/crossfit_folds_K5.rds")
K        <- FD$K
fold_vec <- FD$fold
n_all    <- length(fold_vec)

hdr           <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))
protein_names <- hdr[49:2968]
this_protein  <- protein_names[protein_index]
final_file    <- paste0(OUTDIR, this_protein, "_split_scan.rds")
if (file.exists(final_file)) { cat("already done:", final_file, "\n"); quit(save = "no") }
cat("protein:", this_protein, "\n")

## ---- data, in FD$analysis_eid order (as step 1) ----------------------
df_prot <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE, data.table = FALSE,
                 select = c("eid", "UKBcentre", this_protein))
df_var0 <- fread("data/extracted_pheno_data528_eur.tsv", sep = "\t", header = TRUE, data.table = FALSE)
dis0    <- fread("data/disease_data_level1_eur.tsv", sep = "\t", header = TRUE, data.table = FALSE)
stopifnot(all(dis0$eid == df_var0$eid))
idx <- match(FD$analysis_eid, df_var0$eid); stopifnot(!any(is.na(idx)))
df_var <- df_var0[idx, ]; dis <- dis0[idx, ]; rm(df_var0, dis0); invisible(gc())
df_prot <- df_prot[match(FD$analysis_eid, df_prot$eid), ]
stopifnot(identical(as.character(df_prot$eid), as.character(FD$analysis_eid)),
          identical(as.character(df_var$eid),  as.character(FD$analysis_eid)))
y_all     <- df_prot[, this_protein]
obs_index <- which(!is.na(y_all))

## ---- discovery / validation halves by whole relatedness clusters -----
set.seed(SEED_SPLIT)
blk  <- FD$blocks[sample(length(FD$blocks))]
csz  <- cumsum(lengths(blk))
half <- integer(n_all)
half[unlist(blk[csz <= n_all / 2])] <- 1L
half[unlist(blk[csz >  n_all / 2])] <- 2L
stopifnot(all(half %in% 1:2))
cat(sprintf("discovery n = %d (labeled %d); validation n = %d (labeled %d)\n",
            sum(half == 1), sum(half[obs_index] == 1), sum(half == 2), sum(half[obs_index] == 2)))

## ---- predictors: identical rules to step1_protein_crossfit.R ---------
cand <- colnames(df_var)[8:ncol(df_var)]
cand <- cand[!grepl("^p22009_", cand)]                    # no genetic PCs
X_obs <- as.matrix(df_var[obs_index, cand])
Z_obs <- as.matrix(dis[obs_index, -1])
y_obs <- y_all[obs_index]
fixed_pred <- c("p20116_i0", "p21000_i0", "p31")         # no PRS
sanitize <- function(nm) gsub("[- ]", "_", nm)

screen_fold <- function(train_rows) {
  tp <- match(train_rows, obs_index); stopifnot(!any(is.na(tp)))
  r <- suppressWarnings(cor(X_obs[tp, , drop = FALSE], y_obs[tp], use = "pairwise.complete.obs"))
  r[is.na(r)] <- 0
  cov_pred <- c(colnames(X_obs)[order(abs(r), decreasing = TRUE)[1:100]], fixed_pred)
  yv <- y_obs[tp]
  sig <- vapply(seq_len(ncol(Z_obs)), function(i) {
    zi <- Z_obs[tp, i]; vals <- unique(zi[!is.na(zi)])
    if (length(vals) == 2) -log10(wilcox.test(x = yv[zi == vals[2]], y = yv[zi == vals[1]],
                                              exact = FALSE)$p.value) else 0
  }, numeric(1))
  list(cov_pred = cov_pred, cov_dis = colnames(Z_obs)[order(sig, decreasing = TRUE)[1:10]])
}
build_frame <- function(rows, cov_pred, cov_dis) {
  a <- df_var[rows, cov_pred, drop = FALSE]
  b <- dis[rows, cov_dis, drop = FALSE]
  b[] <- lapply(b, function(x) factor(x, levels = c(0, 1), labels = c("0", "1")))
  out <- cbind(a, b); names(out) <- sanitize(names(out)); out
}

crossfit_half <- function(h) {
  yhat <- rep(NA_real_, n_all)
  for (k in seq_len(K)) {
    train_rows <- intersect(obs_index, which(half == h & fold_vec != k))
    pred_rows  <- which(half == h & fold_vec == k)
    sel <- screen_fold(train_rows)
    tr <- build_frame(train_rows, sel$cov_pred, sel$cov_dis); tr$protein <- y_all[train_rows]
    keep <- vapply(tr[setdiff(names(tr), "protein")], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
    if (any(!keep)) tr <- tr[, c(names(keep)[keep], "protein"), drop = FALSE]
    rf <- ranger(protein ~ ., data = tr, num.trees = 300, save.memory = FALSE,
                 seed = 260611 + 100 * h + k, num.threads = NCORES)
    pr <- build_frame(pred_rows, sel$cov_pred, sel$cov_dis)[, setdiff(names(tr), "protein"), drop = FALSE]
    yhat[pred_rows] <- predict(rf, data = pr)$predictions
    rm(tr, pr, rf); invisible(gc())
  }
  rows <- which(half == h); stopifnot(!any(is.na(yhat[rows])))
  rho <- cor(y_all[intersect(obs_index, rows)], yhat[intersect(obs_index, rows)])
  tick(c("discovery", "validation")[h], ": out-of-fold rho = ", round(rho, 4))
  list(rows = rows, yhat = yhat[rows], rho = rho)
}
CF <- lapply(1:2, crossfit_half)
rm(X_obs, Z_obs); invisible(gc())

## ---- null models per half (same frame construction as step 2) --------
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
cov_id <- c("p21022", paste0("p22009_a", 1:20), "p31", "p22000")

null_half <- function(cf) {
  rows <- cf$rows
  test_df <- df_var[rows, c("eid", cov_id)]
  colnames(test_df) <- c("eid", "age", paste0("geneticPC", 1:20), "sex", "batch")
  test_df$UKBcentre <- df_prot$UKBcentre[rows]
  test_df$age2     <- test_df$age^2
  test_df$age_sex  <- test_df$age * (as.numeric(test_df$sex) - 1)
  test_df$age2_sex <- test_df$age2 * (as.numeric(test_df$sex) - 1)
  test_df$protein     <- y_all[rows]
  test_df$protein_hat <- cf$yhat
  test_df <- INT(test_df, "protein")
  n <- nrow(test_df); r <- rank(test_df$protein_hat)
  test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
  dmm <- model.matrix(~ sex + batch + UKBcentre, data = test_df)
  stopifnot(nrow(dmm) == n)
  tt <- cbind(test_df[, setdiff(colnames(test_df), c("sex", "batch", "UKBcentre"))], dmm)
  tt$`(Intercept)` <- NULL
  colnames(tt) <- gsub("-", "_", colnames(tt))
  GRM <- sgrm[FD$index_all[rows], FD$index_all[rows]]
  mydf <- list(X_all = tt[, setdiff(colnames(tt), c("eid", "protein", "protein_hat", "int", "yhat_int"))],
               S = tt$yhat_int, Y_obs = tt$int, GRM = GRM)
  pS <- SynSurrG_ablation_estimate(mydf)
  vcu <- pS$vc$used
  pS <- add_fastscan_cache(pS, Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(nrow(GRM)))
  pO <- ObsG_ablation_estimate(mydf)
  ## the half's independent subset: FD$independent_indices that fall in the half
  fix_idx <- intersect(which(!is.na(tt$int)), match(intersect(FD$independent_indices, rows), rows))
  list(eid = tt$eid, pS = pS, pO = pO, fix_idx = fix_idx, vc = pS$vc)
}
NM <- lapply(CF, null_half)
tick("null models fitted")
rm(sgrm, df_var, dis); invisible(gc())

## ---- genome-wide scan, both halves per SNP block ----------------------
G0 <- BEDMatrix::BEDMatrix(path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
rs_index <- colnames(G0)
mi <- lapply(NM, function(nm) { m <- match(nm$eid, rownames(G0)); stopifnot(!any(is.na(m))); m })
cn <- as.vector(outer(c("negative_log10_pvalues", "hat_beta", "var_hat_beta"),
                      c("SynSurrG_discovery", "ObsG_discovery", "SynSurrG_validation", "ObsG_validation"),
                      function(a, b) paste0(b, "_", a)))
scan_block <- function(idx) {
  out <- NULL
  for (h in 1:2) {
    Gm <- as.matrix(G0[mi[[h]], idx]); Gm <- replace(Gm, is.na(Gm), 0)
    Gm <- fix_constant_columns(Gm, NM[[h]]$fix_idx)
    r1 <- score_test_SynSurrG_fast(Gm, NM[[h]]$pS, label = "SynSurrG")
    r2 <- score_test_ObsG_multiply(Gm, NM[[h]]$pO)
    out <- cbind(out, r1$negative_log10_pval_SynSurrG, r1$hat_beta_SynSurrG, r1$var_hat_beta_SynSurrG,
                      r2$negative_log10_pval_ObsG,     r2$hat_beta_ObsG,     r2$var_hat_beta_ObsG)
  }
  dimnames(out) <- list(colnames(G0)[idx], cn); out
}
seg_list <- split(seq_along(rs_index), cut(seq_along(rs_index), breaks = NSEG, labels = FALSE))
seg_file <- function(s) paste0(OUTDIR, this_protein, "_seg", s, "of", NSEG, RUN_TAG, ".rds")
for (s in seq_len(NSEG)) {
  if (file.exists(seg_file(s))) { tick("segment ", s, " already done"); next }
  blocks <- split(seg_list[[s]], ceiling(seq_along(seg_list[[s]]) / 200))
  res <- parallel::mclapply(blocks, function(i) tryCatch(scan_block(i), error = function(e) e),
                            mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) !is.matrix(x), logical(1))
  if (any(bad)) stop("segment ", s, ": ", sum(bad), " block(s) failed; first: ", conditionMessage(res[[which(bad)[1]]]))
  saveRDS(do.call(rbind, res), seg_file(s))
  tick("segment ", s, "/", NSEG, " saved")
}
results <- do.call(rbind, lapply(seq_len(NSEG), function(s) readRDS(seg_file(s))))
stopifnot(identical(rownames(results), rs_index))
attr(results, "rho")   <- c(discovery = CF[[1]]$rho, validation = CF[[2]]$rho)
attr(results, "n")     <- c(discovery = length(CF[[1]]$rows), validation = length(CF[[2]]$rows))
attr(results, "n_obs") <- c(discovery = sum(!is.na(y_all[CF[[1]]$rows])), validation = sum(!is.na(y_all[CF[[2]]$rows])))
attr(results, "vc")    <- list(discovery = NM[[1]]$vc, validation = NM[[2]]$vc)
attr(results, "seed_split") <- SEED_SPLIT
attr(results, "code")  <- "split_validation_crossfit.R (2026-10-13)"
saveRDS(results, final_file)
invisible(file.remove(vapply(seq_len(NSEG), seg_file, "")))
tick("done: ", final_file)

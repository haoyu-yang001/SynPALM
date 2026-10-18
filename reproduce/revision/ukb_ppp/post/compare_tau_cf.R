## Supplementary figure "twoSGWAS" (paper/figs/SynGWASalphaG.pdf), revision:
## Syn-PALM -log10 P with the variance components (VCs) estimated once under
## the global null (beta_G = 0, alpha_G = 0; production setting) versus VCs
## re-estimated for each SNP under beta_G = 0 only (the SNP is a covariate of
## the synthetic phenotype when the VCs are estimated), for the 2,000 SNPs most
## strongly associated with the synthetic phenotype of one well-predicted protein
## (default LEP, out-of-fold rho 0.884; the submitted figure used FOLR3 with
## PRS-based Yhat, which is no longer predictable without PRS: rho 0.026).
## Replaces code/compare_taus.R. Data and null model exactly as
## code/synpalm_main_code/step2_protein_fast.R (step-1 cache
## <p>_synthetic_phenotype_cf_noPRS_noPC.rds, crossfit_folds_K5.rds, sparse GRM,
## covariates, INT(Y), rank-INT(Yhat), v2 array genotypes with NA -> 0 and
## fix_constant_columns, SynSurrG_ablation_estimate + add_fastscan_cache +
## score_test_SynSurrG_fast).
## SNP-specific VCs: as "SynSurrG_each" in code/compare_taus.R, i.e.
## SynSurrG_eachG_step1(mydf, G) of code/funcs_v2.R, which is
## SynSurrG_ablation_estimate with the SNP added ONLY to the regression of Yhat
## on the covariates whose residuals feed the method-of-moments estimates
## (tau_S2, sigma_S2, tau_TS, sigma_TS); the regression of Y keeps beta_G = 0, and
## the mean model of the score test (X_all) does not contain the SNP, which is
## then scored with the refitted parameters (alpha_G is estimated inside the
## test as in production). The covariate copy of the SNP is mean-imputed.
## Selection of the 2,000 SNPs: linear regression of rank-INT(Yhat) on the
## same covariates (with intercept) in all n analysed individuals, one SNP at a
## time (NA -> 0, as the big_univLinReg call that made the old
## folr3_all_hatYlog10lpvals.rds, commented out in code/compare_taus.R),
## autosomes only, MHC chr6:25.5-34.0 Mb (GRCh38) excluded, top 2,000 by P.
##
##   Rscript code/post_cf3/compare_tau_cf.R yhatgwas <array> <n_array>  # GWAS of Yhat, one SNP segment
##   Rscript code/post_cf3/compare_tau_cf.R top <n_array>               # merge, pick top 2,000
##   Rscript code/post_cf3/compare_tau_cf.R refit <array> <n_array>     # both P for one share of the 2,000
## (SLURM driver: code/post_cf3/submit_compare_tau_cf.sh;
##  merge + figure: code/post_cf3/compare_tau_cf_plot.R)
## env: SYNPALM_CT_PROTEIN protein index (default 1571 = lep)
##      SYNPALM_CT_NTOP    number of SNPs (default 2000)
## Outputs in results/compare_tau_cf_noPRS_noPC/:
##   temp/<p>_yhatgwas_seg<a>of<N>.rds   per-segment Yhat GWAS (beta, se, -log10 P)
##   <p>_yhatgwas_cf_noPRS_noPC.rds      merged Yhat GWAS (autosomes, no MHC)
##   <p>_top2000_yhat_snps.rds           selected SNPs, in order of Yhat -log10 P
##   temp/<p>_refit_arr<a>of<N>.rds      per-array: global-null and SNP-specific
##                                       -log10 P, beta, var(beta), refitted VCs
## A finished segment / array (and, within an array, each finished SNP via
## temp/<p>_refit_arr<a>of<N>_partial.rds) is skipped on restart.
## =====================================================================

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(Matrix) })
setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
S1_TAG  <- "_cf_noPRS_noPC"
t_start <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")

args  <- commandArgs(trailingOnly = TRUE)
MODE  <- args[1]
stopifnot(MODE %in% c("yhatgwas", "top", "refit"))
ARR   <- if (length(args) >= 2) as.integer(args[2]) else 1L
NARR  <- if (length(args) >= 3) as.integer(args[3]) else 1L
protein_index <- as.integer(Sys.getenv("SYNPALM_CT_PROTEIN", "1571"))
NTOP  <- as.integer(Sys.getenv("SYNPALM_CT_NTOP", "2000"))
BED   <- "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed"
OUT   <- "results/compare_tau_cf_noPRS_noPC/"
TMP   <- paste0(OUT, "temp/")
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

protein_names <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
this_protein  <- protein_names[protein_index]
gwas_file <- paste0(OUT, this_protein, "_yhatgwas", S1_TAG, ".rds")
top_file  <- paste0(OUT, this_protein, "_top", NTOP, "_yhat_snps.rds")

## autosomal, non-MHC SNPs of the v2 array (same annotation as the paper figures)
snp_set <- function() {
  ann <- readRDS("data/snp_annot_v2geno.rds")
  ann <- ann[chr %in% as.character(1:22) & !(chr == "6" & POS38 >= 25.5e6 & POS38 <= 34e6)]
  ann
}

## ---------------------------------------------------------------- top
if (MODE == "top") {
  ## NARR here = number of yhatgwas segments
  NSEG <- if (length(args) >= 2) as.integer(args[2]) else 20L
  fs <- sprintf("%s%s_yhatgwas_seg%dof%d.rds", TMP, this_protein, seq_len(NSEG), NSEG)
  miss <- fs[!file.exists(fs)]
  if (length(miss)) stop(length(miss), " Yhat-GWAS segment(s) missing, first: ", miss[1])
  gw <- do.call(rbind, lapply(fs, readRDS))
  ann <- snp_set()
  stopifnot(!anyDuplicated(rownames(gw)), all(rownames(gw) %in% ann$rsid))
  gw <- data.table(rsid = rownames(gw), as.data.table(gw))
  gw <- merge(ann[, .(rsid, chr, POS38)], gw, by = "rsid", sort = FALSE)
  attr(gw, "method") <- "lm of rank-INT(Yhat) on covariates (+ intercept), all individuals, NA->0"
  saveRDS(gw, gwas_file)
  setorder(gw, -yhat_log10p)
  top <- gw[seq_len(NTOP)]
  saveRDS(top, top_file)
  tick(nrow(gw), " SNPs; top ", NTOP, " selected: -log10 P(Yhat) from ",
       sprintf("%.1f to %.1f", min(top$yhat_log10p), max(top$yhat_log10p)),
       "; chromosomes: ", paste(names(table(top$chr)), table(top$chr), sep = ":", collapse = " "))
  for (s in fs) file.remove(s)
  quit(save = "no")
}

## ---------------------------------------------------------------- data (as step2_protein_fast.R)
seg_out <- sprintf("%s%s_yhatgwas_seg%dof%d.rds", TMP, this_protein, ARR, NARR)
arr_out <- sprintf("%s%s_refit_arr%dof%d.rds", TMP, this_protein, ARR, NARR)
arr_tmp <- sub("\\.rds$", "_partial.rds", arr_out)
if (MODE == "yhatgwas" && file.exists(seg_out)) { tick("segment already done"); quit(save = "no") }
if (MODE == "refit"    && (file.exists(arr_out) ||
    file.exists(paste0(OUT, this_protein, "_TOTAL_compare_tau_cf_noPRS_noPC.rds")))) {
  tick("array already done (or already merged by compare_tau_cf_plot.R)"); quit(save = "no") }

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
cat(sprintf("protein %s (index %d): rho_oof = %.4f, n_obs = %d\n",
            this_protein, protein_index, S1$rho_oof, S1$n_obs))

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

G0 <- BEDMatrix::BEDMatrix(path = BED, simple_names = TRUE)
matched_indices <- match(test_df_transformed$eid, rownames(G0))
if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
fix_idx <- intersect(obs_protein_index, FD$independent_indices)
tick("data ready (n = ", n, ", n_obs = ", length(obs_protein_index), ")")

## ---------------------------------------------------------------- yhatgwas
if (MODE == "yhatgwas") {
  ann  <- snp_set()
  snps <- intersect(colnames(G0), ann$rsid)          # BED order
  segs <- split(snps, cut(seq_along(snps), breaks = NARR, labels = FALSE))
  mine <- segs[[ARR]]
  X  <- cbind(1, as.matrix(test_df_transformed[, X_all_names]))
  qx <- qr(X); Q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
  s  <- test_df_transformed$yhat_int
  rs <- as.numeric(s - Q %*% crossprod(Q, s)); rss0 <- sum(rs^2)
  df_res <- n - qx$rank - 1
  set.seed(1000 + ARR)                                # fix_constant_columns samples rows
  blocks <- split(mine, ceiling(seq_along(mine) / 1000))
  res <- lapply(blocks, function(b) {
    Gm <- as.matrix(G0[matched_indices, b, drop = FALSE])
    Gm <- replace(Gm, is.na(Gm), 0)
    Gm <- fix_constant_columns(Gm, fix_idx)
    Rg <- Gm - Q %*% crossprod(Q, Gm)
    gg <- colSums(Rg^2); gs <- as.numeric(crossprod(Rg, rs))
    beta <- gs / gg
    se   <- sqrt((rss0 - beta * gs) / df_res / gg)
    lp   <- -(pt(abs(beta / se), df_res, lower.tail = FALSE, log.p = TRUE) + log(2)) / log(10)
    cbind(yhat_beta = beta, yhat_se = se, yhat_log10p = lp)
  })
  out <- do.call(rbind, res); rownames(out) <- mine
  saveRDS(out, seg_out)
  tick("segment ", ARR, "/", NARR, " saved (", nrow(out), " SNPs)")
  quit(save = "no")
}

## ---------------------------------------------------------------- refit
top  <- readRDS(top_file)
alls <- top$rsid
mine <- split(alls, cut(seq_along(alls), breaks = NARR, labels = FALSE))[[ARR]]

sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
GRM <- sgrm[FD$index_all, FD$index_all]
rm(sgrm); invisible(gc())
mydf <- list(X_all = test_df_transformed[, X_all_names],
             S     = test_df_transformed$yhat_int,
             Y_obs = test_df_transformed$int,
             GRM   = GRM)

set.seed(2000 + ARR)                                  # fix_constant_columns samples rows
Gm  <- as.matrix(G0[matched_indices, mine, drop = FALSE])
Gmi <- apply(Gm, 2, function(g) { g[is.na(g)] <- mean(g, na.rm = TRUE); g })   # covariate copy
Gm  <- replace(Gm, is.na(Gm), 0)
Gm  <- fix_constant_columns(Gm, fix_idx)               # scored copy, as production
tick(length(mine), " SNPs of array ", ARR, "/", NARR, " loaded")

vc_row <- function(vc) {        # used (PSD-constrained) and raw method-of-moments VCs
  nm <- c("tau_T2", "tau_TS", "tau_S2", "sigma_T2", "sigma_TS", "sigma_S2")
  c(setNames(vc$used[nm], nm), setNames(vc$raw[nm], paste0(nm, "_raw")),
    projected_G = as.numeric(isTRUE(vc$projected_G)), projected_E = as.numeric(isTRUE(vc$projected_E)),
    rho_g = as.numeric(vc$rho_g), rho_e = as.numeric(vc$rho_e))
}
cache <- function(p) {
  vcu <- p$vc$used
  add_fastscan_cache(p, Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(nrow(GRM)))
}

## global null: once per array (all 2,000 SNPs share it)
P0 <- cache(SynSurrG_ablation_estimate(mydf))
r0 <- score_test_SynSurrG_fast(Gm, P0, label = "SynSurrG")
vc0 <- vc_row(P0$vc); rm(P0); invisible(gc())
tick("global-null model + scan")

## SNP-specific refits; with SLURM_CPUS_PER_TASK > 1 they run in batches on
## forked workers (the deterministic refits do not depend on the batching)
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
refit_one <- function(j) {
  t0 <- Sys.time()
  Pj <- cache(SynSurrG_eachG_step1(mydf, G = Gmi[, j]))
  rj <- score_test_SynSurrG_fast(Gm[, j, drop = FALSE], Pj, label = "SynSurrG")
  c(SynSurrG_each_negative_log10_pvalues = rj$negative_log10_pval_SynSurrG,
    SynSurrG_each_hat_beta = rj$hat_beta_SynSurrG,
    SynSurrG_each_var_hat_beta = rj$var_hat_beta_SynSurrG,
    setNames(vc_row(Pj$vc), paste0("each_", names(vc0))),
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
done <- if (file.exists(arr_tmp)) readRDS(arr_tmp) else list()
todo <- which(!mine %in% names(done))
for (b in split(todo, ceiling(seq_along(todo) / NCORES))) {
  rb <- parallel::mclapply(b, function(j) tryCatch(refit_one(j), error = function(e) e),
                           mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(rb, function(x) inherits(x, "error") || is.null(x), logical(1))
  if (any(bad)) stop("refit failed for ", paste(mine[b][bad], collapse = ","), ": ",
                     if (inherits(rb[[which(bad)[1]]], "error")) conditionMessage(rb[[which(bad)[1]]]) else "NULL (worker killed?)")
  for (k in seq_along(b)) done[[mine[b[k]]]] <- rb[[k]]
  saveRDS(done, arr_tmp)
  for (j in b) tick("SNP ", j, "/", length(mine), " ", mine[j], ": -log10P global ",
                    sprintf("%.3f", r0$negative_log10_pval_SynSurrG[j]), ", SNP-specific ",
                    sprintf("%.3f", done[[mine[j]]][1]))
  invisible(gc())
}
each <- do.call(rbind, done[mine])
res  <- cbind(SynSurrG_negative_log10_pvalues = r0$negative_log10_pval_SynSurrG,
              SynSurrG_hat_beta = r0$hat_beta_SynSurrG,
              SynSurrG_var_hat_beta = r0$var_hat_beta_SynSurrG,
              each)
rownames(res) <- mine
attr(res, "vc_global") <- vc0
saveRDS(res, arr_out)
invisible(file.remove(arr_tmp))
tick("array ", ARR, "/", NARR, " saved")

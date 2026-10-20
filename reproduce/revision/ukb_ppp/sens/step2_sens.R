## SENSITIVITY COPY (revision 2026-10-19; reply to Reviewer 1, Comments 4 and 8) of
## code/synpalm_main_code/step2_protein_fast.R: identical per-SNP computation, except that
##   SENS_SNP_FILE  only the rsids listed in this file are scanned (one per line)
##   only Syn-PALM and Obs-GWAS are computed (SynSurr / Obs columns are NA)
## and the output is <protein><S1_TAG>_scan_subset.rds
## =====================================================================
## STEP 2 of 2 (production) -- genome-wide scan of ALL SNPs for ONE protein.
## Null models are fitted once; the genome is cut into NSEG segments, each
## segment is scanned in blocks spread over SLURM_CPUS_PER_TASK forked workers
## and saved on its own, so a rerun skips finished segments.
## Same per-SNP computation as step2_scan_fast.R (validated against funcs_v2).
##
##   Rscript step2_protein_fast.R <protein_index> [NSEG]
##
## (header of step2_scan_fast.R follows)
## Reads the cross-fitted synthetic phenotype cached by step 1 and scans with
## the vectorised score test in code/fastscan/funcs_fastscan.R (same output as
## score_test_SynSurrG_multiply, without the per-SNP loop).
##
##   Rscript step2_scan_fast.R <protein_index> <chunk_array> <arraynum>
##
## env: SYNPALM_S1_TAG   step-1 tag to read (default _cf_noPRS_noPC; a
##                       permutation run is e.g. _cf_noPRS_noPC_PERMUTED1)
##      SYNPALM_COND_SNPS comma-separated rsids added as covariates to every null
##                       model (conditional scan, e.g. a cis-pQTL that dominates the
##                       protein); their own rows are copied from SYNPALM_UNCOND_FILE,
##                       the unconditional scan of the same protein
##      SYNPALM_VALIDATE =1 also runs the funcs_v2 per-SNP functions on the
##                       first 200 SNPs and prints the largest differences
## =====================================================================

library(data.table)
library(dplyr)
library(Matrix)

setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")
source("code/fastscan/funcs_fastscan.R")
S1_TAG   <- Sys.getenv("SYNPALM_S1_TAG", "_cf_noPRS_noPC")
VALIDATE <- identical(Sys.getenv("SYNPALM_VALIDATE", "0"), "1")
t_start  <- Sys.time()
tick <- function(...) cat(sprintf("[%6.1f s] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")

args          <- commandArgs(trailingOnly = TRUE)
protein_index <- as.integer(args[1])
NSEG          <- if (length(args) >= 2) as.integer(args[2]) else 10L
NCORES        <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])

## covariates only -- the 2,920-column protein table is not needed here
df_allUKB_eur <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE,
                       data.table = FALSE, select = c("eid", "UKBcentre"))
df_all_variables_eur_orignal <- fread("data/extracted_pheno_data528_eur.tsv",
                                      sep = "\t", header = TRUE, data.table = FALSE)
protein_names <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
this_protein  <- protein_names[protein_index]

FD <- readRDS("data/crossfit_folds_K5.rds")
idx <- match(FD$analysis_eid, df_all_variables_eur_orignal$eid)
df_all_variables_eur <- df_all_variables_eur_orignal[idx, ]
rm(df_all_variables_eur_orignal); gc()

outdir  <- paste0("results/allukb_protein_pQTL_european_discovery/", this_protein, "/")
s1_file <- paste0(outdir, this_protein, "_synthetic_phenotype", S1_TAG, ".rds")
if (!file.exists(s1_file))
  stop("step-1 output missing for ", this_protein, ": ", s1_file,
       "\nRun step1_protein_crossfit.R for this protein first.")
S1 <- readRDS(s1_file)

## Guards. A permuted synthetic phenotype must never reach the scan, and the
## PRS setting used in step 1 is recorded here so the log says which it was.
stopifnot(identical(as.character(S1$eid), as.character(FD$analysis_eid)))
if (isTRUE(S1$permuted) != grepl("_PERMUTED", S1_TAG))
  stop("step-1 cache ", s1_file, " has permuted = ", S1$permuted, " but S1_TAG = ", S1_TAG)
if (!identical(S1$use_prs, FALSE) || !identical(S1$use_pcs, FALSE))
  stop("step-1 cache ", s1_file, " was built with PRS or PCs in the forest")
cat(sprintf("using step-1 cache %s: rho_oof = %.4f, n_obs = %d, use_prs = %s, use_pcs = %s\n",
            basename(s1_file), S1$rho_oof, S1$n_obs, S1$use_prs, S1$use_pcs))

## the observed protein for this protein only
y_all <- fread("data/df_allUKB_eur1.tsv", sep = "\t", header = TRUE, data.table = FALSE,
               select = c("eid", this_protein))
y_all <- y_all[match(FD$analysis_eid, y_all$eid), this_protein]

GRM                 <- sgrm[FD$index_all, FD$index_all]
independent_indices <- FD$independent_indices

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

COND_SNPS <- Filter(nzchar, strsplit(Sys.getenv("SYNPALM_COND_SNPS", ""), ",")[[1]])
if (length(COND_SNPS)) {
  Gc <- BEDMatrix::BEDMatrix(path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed", simple_names = TRUE)
  gc_mat <- as.matrix(Gc[match(test_df_transformed$eid, rownames(Gc)), COND_SNPS, drop = FALSE])
  for (j in seq_along(COND_SNPS)) {
    g <- gc_mat[, j]; g[is.na(g)] <- mean(g, na.rm = TRUE)
    test_df_transformed[[paste0("cond_", COND_SNPS[j])]] <- g
  }
  rm(Gc, gc_mat)
  tick("conditioning on ", paste(COND_SNPS, collapse = ", "))
}
obs_protein_index <- which(!is.na(test_df_transformed$protein))
X_all_names <- setdiff(colnames(test_df_transformed),
                       c("eid", "protein", "protein_hat", "int", "yhat_int"))
mydf <- list(X_all = test_df_transformed[, X_all_names],
             S     = test_df_transformed$yhat_int,
             Y_obs = test_df_transformed$int,
             GRM   = GRM)

tick("data ready")
SynSurrG_step1_pars <- SynSurrG_ablation_estimate(mydf)
vcu <- SynSurrG_step1_pars$vc$used
SynSurrG_step1_pars <- add_fastscan_cache(SynSurrG_step1_pars,
                         Sigma22 = vcu[["tau_S2"]] * GRM + vcu[["sigma_S2"]] * Diagonal(nrow(GRM)))
tick("Syn-PALM null model")
ObsG_step1_pars     <- ObsG_ablation_estimate(mydf)
tick("other null models")
## the forked workers only need the parameter lists: drop the big inputs first
rm(sgrm, GRM, mydf); invisible(gc())

G0 <- BEDMatrix::BEDMatrix(path = "genotype_data/ukb_cal_allchr_v2_whiteb_QC.bed",
                           simple_names = TRUE)
geno_samp.id  <- rownames(G0)
rs_index      <- intersect(colnames(G0), readLines(Sys.getenv("SENS_SNP_FILE")))
tick("scanning ", length(rs_index), " SNPs of ", Sys.getenv("SENS_SNP_FILE"))
seg_list      <- if (NSEG == 1) list(rs_index) else split(rs_index, cut(seq_along(rs_index), breaks = NSEG, labels = FALSE))
chunk_size    <- as.integer(Sys.getenv("SYNPALM_CHUNK", "200"))
matched_indices <- match(test_df_transformed$eid, geno_samp.id)
if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
fix_idx <- intersect(obs_protein_index, independent_indices)
seg_file <- function(s) paste0(outdir, this_protein, "_seg", s, "of", NSEG, S1_TAG, "_scan.rds")
final_file <- paste0(outdir, this_protein, S1_TAG, "_scan_subset.rds")
cn <- c("SynSurrG_negative_log10_pvalues","SynSurrG_hat_beta","SynSurrG_var_hat_beta",
        "ObsG_negative_log10_pvalues","ObsG_hat_beta","ObsG_var_hat_beta",
        "SynSurr_negative_log10_pvalues","SynSurr_hat_beta","SynSurr_var_hat_beta",
        "Obs_negative_log10_pvalues","Obs_hat_beta","Obs_var_hat_beta")

scan_block <- function(idx) {
  Gm <- as.matrix(G0[matched_indices, idx])
  Gm <- replace(Gm, is.na(Gm), 0)
  Gm[, colnames(Gm) %in% COND_SNPS] <- 0     # tested later from the unconditional scan
  Gm <- fix_constant_columns(Gm, fix_idx)
  r1 <- score_test_SynSurrG_fast(Gm, SynSurrG_step1_pars, label = "SynSurrG")
  r2 <- score_test_ObsG_multiply(Gm, ObsG_step1_pars)
  na <- rep(NA_real_, ncol(Gm))
  r3 <- list(negative_log10_pval_SynSurr = na, hat_beta_SynSurr = na, var_hat_beta_SynSurr = na)
  r4 <- list(negative_log10_pval_Obs = na, hat_beta_Obs = na, var_hat_beta_Obs = na)
  out <- cbind(r1$negative_log10_pval_SynSurrG, r1$hat_beta_SynSurrG, r1$var_hat_beta_SynSurrG,
               r2$negative_log10_pval_ObsG,     r2$hat_beta_ObsG,     r2$var_hat_beta_ObsG,
               r3$negative_log10_pval_SynSurr,  r3$hat_beta_SynSurr,  r3$var_hat_beta_SynSurr,
               r4$negative_log10_pval_Obs,      r4$hat_beta_Obs,      r4$var_hat_beta_Obs)
  dimnames(out) <- list(colnames(Gm), cn)
  out
}

if (!file.exists(final_file)) {
  for (s in seq_len(NSEG)) {
    if (file.exists(seg_file(s))) { tick("segment ", s, " already done"); next }
    snps   <- seg_list[[s]]
    blocks <- split(snps, ceiling(seq_along(snps) / chunk_size))
    res <- parallel::mclapply(blocks, function(idx) tryCatch(scan_block(idx), error = function(e) e),
                              mc.cores = NCORES, mc.preschedule = FALSE)
    bad <- vapply(res, function(x) !is.matrix(x), logical(1))
    if (any(bad)) stop("segment ", s, ": ", sum(bad), " block(s) failed; first: ",
                       conditionMessage(res[[which(bad)[1]]]))
    seg <- do.call(rbind, res)
    stopifnot(identical(rownames(seg), snps))
    saveRDS(seg, seg_file(s))
    tick("segment ", s, "/", NSEG, " saved (", nrow(seg), " SNPs)")
  }
  results <- do.call(rbind, lapply(seq_len(NSEG), function(s) readRDS(seg_file(s))))
  stopifnot(identical(rownames(results), rs_index))
  attr(results, "rho_oof")  <- S1$rho_oof
  attr(results, "n_obs")    <- S1$n_obs
  attr(results, "use_prs")  <- S1$use_prs
  attr(results, "use_pcs")  <- S1$use_pcs
  attr(results, "s1_tag")   <- S1_TAG
  attr(results, "n_total")  <- nrow(test_df_transformed)
  attr(results, "vc")       <- SynSurrG_step1_pars$vc
  if (length(COND_SNPS)) {
    unc <- readRDS(Sys.getenv("SYNPALM_UNCOND_FILE"))
    results[COND_SNPS, ] <- unc[COND_SNPS, colnames(results)]
    attr(results, "cond_snps") <- COND_SNPS
    attr(results, "vc_uncond") <- attr(unc, "vc")
  }
  attr(results, "code")     <- "step2_protein_fast.R + funcs_fastscan.R (2026-10-13)"
  saveRDS(results, final_file)
  invisible(file.remove(vapply(seq_len(NSEG), seg_file, "")))
}
tick("Protein ", this_protein, " done: ", final_file)

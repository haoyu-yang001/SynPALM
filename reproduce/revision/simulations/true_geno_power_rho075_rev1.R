## Rerun of the rho = 0.75 settings of the true-genotype power simulation (Supp
## Figs power_Syn-PALM_ObsG_true_geno, power_ratio_Syn-PALM_vs_{obsG,SynSurr}_true_geno;
## simulation_true_geno_power_chr_v2.R), whose genetic covariance broke
## Cauchy-Schwarz (Reviewer 2, comment 5). Same data generation (UKB sparse GRM,
## first n eligible White-British eids, DGP_GRMr_nullmodel_step1-2 and
## DGP_GRMr_gmodel_step3 with the real genotype of the SNP), same SNP sample
## (with_seed(241150) sample of 20,000 genotyped SNPs, split into 1000 arrays of 20),
## same per-SNP phenotype and per-SNP null-model refit for all six methods, same
## *_single tests and output columns; changes:
##   - rho = 0.75 only, with sim_cov_TS() (genetic correlation 0.99, total rho 0.75);
##   - autosomal SNPs only: the non-autosomal SNPs of the 20,000 sample (chromosome
##     X, Y, XY, MT in the .bim) are skipped, so rows line up with the autosomal
##     rows of the submitted rho < 0.75 results (same SNPs, same order);
##   - every random draw is seeded: one seed per (array, missing rate, pve_g),
##     saved in attr(, "settings"); the fix_constant_columns draw per (array, miss);
##   - the genotype block of an array is read once per missing rate (the submitted
##     script re-read the same block for every pve_g);
##   - one output file per (missing rate, pve_g, array), skipped if it exists, so a
##     preempted / requeued task resumes where it stopped;
##   - from 2026-10-14 only the three methods shown in the figures are computed
##     (Syn-PALM = SynSurrG, Obs-GWAS = ObsG, SynSurr); the OracleG, Oracle and Obs
##     columns are NA (these tests draw no random numbers, so the other columns are
##     unchanged; files written earlier have all 18 columns; attr(, "methods") says which).
##   Rscript true_geno_power_rho075_rev1.R <array> [miss list, e.g. 0,0.5] [max_snps] [out_dir] [rho=0.75]
## Revision 2 (2026-10-16): rho = 0, 0.25, 0.5 (kept from the submission) are rerun with
## this script too (argument rho); output in results_rev1/true_geno_power_rho<100 rho>/,
## phenotype seeds offset by 1e7 x (1, 2, 3) for rho = 0, 0.25, 0.5; the
## fix_constant_columns seed does not depend on rho, so all rho use the same genotypes.
## Output: simulation_code/results_rev1/true_geno_power_rho075/
##   missing_<m>_rho_0.75_her_<h>_chunkarray_<a>_negative_log10_pvals_for_power_rev1.rds
##   (18 columns as the submitted files; rsids in attr(, "snps")).
setwd("<PROJECT_DIR>")
suppressMessages(source("code/funcs_v2.R"))
options(synsurrg.verbose = FALSE)
library(withr)

arraynum <- 1000L                                   # as submitted: 20,000 SNPs / 1000 arrays
args     <- commandArgs(trailingOnly = TRUE)
array_id <- as.integer(args[1])
miss_all <- c(0, 0.25, 0.5, 0.75, 0.9)
miss_vec <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(strsplit(args[2], ",")[[1]]) else miss_all
max_snps <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else NA_integer_
rho      <- if (length(args) >= 5 && nzchar(args[5])) as.numeric(args[5]) else 0.75
rho_off  <- c("0.75" = 0L, "0" = 10000000L, "0.25" = 20000000L, "0.5" = 30000000L)[[as.character(rho)]]
OUT      <- if (length(args) >= 4 && nzchar(args[4])) args[4] else sprintf("simulation_code/results_rev1/true_geno_power_rho%03d/", round(100 * rho))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

snp_heritability_vec <- seq(0.0001, 0.0007, 0.0001)
n_obs <- 3e4
tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5
pve_x <- 0
cols <- c("OracleG_negative_log10_pval","OracleG_hat_beta","OracleG_var_hat_beta",
          "ObsG_negative_log10_pval","ObsG_hat_beta","ObsG_var_hat_beta",
          "SynSurrG_negative_log10_pval","SynSurrG_hat_beta","SynSurrG_var_hat_beta",
          "Oracle_negative_log10_pval","Oracle_hat_beta","Oracle_var_hat_beta",
          "Obs_negative_log10_pval","Obs_hat_beta","Obs_var_hat_beta",
          "SynSurr_negative_log10_pval","SynSurr_hat_beta","SynSurr_var_hat_beta")

out_file <- function(miss, pve_g) sprintf("%smissing_%s_rho_%s_her_%s_chunkarray_%d_negative_log10_pvals_for_power_rev1.rds",
                                          OUT, miss, rho, pve_g, array_id)
grid <- expand.grid(pve_g = snp_heritability_vec, miss = miss_vec)
grid <- grid[!file.exists(mapply(out_file, grid$miss, grid$pve_g)), ]
if (!nrow(grid)) { cat("all outputs exist for array", array_id, "\n"); quit(save = "no") }

t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

### GRM, sample and genotypes: as in the submitted script ###################
sgrm0 <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
rm(sGRM)
sgrm0 <- sgrm0*2
colnames(sgrm0) <- sapply(colnames(sgrm0), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm0) <- sapply(rownames(sgrm0), function(x) strsplit(x, "_")[[1]][1])
blocks <- find_blocks_vectorized(sgrm0)
block_num <- sapply(blocks, length)
small_members <- unlist(blocks[which(block_num < 500)])
sgrm_bf <- sgrm0[small_members, small_members]
hh1 <- find_blocks_vectorized(sgrm_bf)
perm_idx <- with_seed(255155, { sample(length(hh1)) })
idx <- unlist(hh1[perm_idx], use.names = FALSE)
sgrm <- sgrm_bf[idx, idx, drop = FALSE]
rm(sgrm0, sgrm_bf, blocks, hh1); invisible(gc())
eids_sGRM <- colnames(sgrm)

direction <- "<UKB_DIR>/Pharma_Proteomics_Project_2923_proteins_data/"
df_phenos <- fread(paste0(direction, "/phenotypes/extracted_pheno_data.tsv"), data.table = F)
colnames(df_phenos) <- gsub("^participant\\.", "", colnames(df_phenos))
df_phenos_colnames <- fread(paste0(direction, "/phenotypes/extracted_pheno_data_meaning.csv"), data.table = F)
colnames(df_phenos) <- df_phenos_colnames[-1, 2]
subeids0 <- df_phenos$eid[which(is.na(df_phenos$`UKB-PPP Consortium selected participant | Instance 0`) &
                                  is.na(df_phenos$`Case-control status for COVID19 imaging repeat | Instance 3 | Array 0`) &
                                  !is.na(df_phenos$`Genetic ethnic grouping`))]
subeids <- intersect(eids_sGRM, subeids0)
rm(df_phenos); invisible(gc())

geno_prefix <- "genotype_data/ukb_cal_allchr_v2_whiteb_QC"
G0 <- BEDMatrix::BEDMatrix(path = paste0(geno_prefix, ".bed"), simple_names = T)
geno_samp.id <- rownames(G0)
bim <- fread(paste0(geno_prefix, ".bim"), header = FALSE, data.table = FALSE)
stopifnot(identical(colnames(G0), bim$V2))
rs_index <- colnames(G0)
select_index <- with_seed(241150, sample(1:length(rs_index), 200*100, replace = F))   # submitted SNP sample
rs_index <- rs_index[select_index]
rs_index_list <- split(rs_index, cut(seq_along(rs_index), breaks = arraynum, labels = FALSE))
is_auto <- setNames(bim$V1 %in% as.character(1:22), bim$V2)
rs_array <- rs_index_list[[array_id]]
rs_array <- rs_array[is_auto[rs_array]]                                                 # autosomes only
tick("array ", array_id, ": ", length(rs_array), " of ", length(rs_index_list[[array_id]]), " sampled SNPs autosomal")
if (!is.na(max_snps)) rs_array <- head(rs_array, max_snps)

for (miss in unique(grid$miss)) {
  mi <- match(miss, miss_all)
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  index <- match(subeids[1:n_all], eids_sGRM)
  GRM <- sgrm[index, index]
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all

  hh <- find_blocks_vectorized(GRM)
  independent_indices <- with_seed(25324, {
    sapply(hh, function(x) x[sample.int(length(x), 1)])
  })

  pre_vars_step1 <- DGP_GRMr_nullmodel_step1(n_obs, miss, tauT, tauS, sigmaT, sigmaS, GRM)
  pre_vars_step2 <- DGP_GRMr_nullmodel_step2(n_obs, miss, tauT, tauS, sigmaT, sigmaS, rho, GRM, pre_vars_step1)
  rm(pre_vars_step1); invisible(gc())

  matched_indices <- match(subeids[1:n_all], geno_samp.id)
  if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
  seed_g <- 20261114L + 1000L * array_id + 10L * mi
  Genotype_matrix <- as.matrix(G0[matched_indices, rs_array, drop = FALSE])
  Genotype_matrix <- replace(Genotype_matrix, is.na(Genotype_matrix), 0)
  Genotype_matrix <- with_seed(seed_g, fix_constant_columns(Genotype_matrix, intersect(1:n_obs, independent_indices)))
  tick("miss ", miss, ": DGP prepared (n = ", n_all, ")")

  for (pve_g in grid$pve_g[grid$miss == miss]) {
    hi <- match(pve_g, snp_heritability_vec)
    seed <- seed_g + hi + rho_off
    set.seed(seed)
    res <- matrix(NA_real_, length(rs_array), length(cols), dimnames = list(NULL, cols))
    for (ss in seq_along(rs_array)) {
      mydf <- DGP_GRMr_gmodel_step3(n_obs, miss, tauT = tauT, sigmaT = sigmaT, pve_x = pve_x, pve_g = pve_g,
                                    maf = NA, kinsetting = NA, Genotype_matrix[, ss], pre_vars_step2, GRM)
      G_all <- mydf$G_all
      na3 <- function(m) setNames(list(NA_real_, NA_real_, NA_real_),
                                  paste0(c("negative_log10_pval_", "hat_beta_", "var_hat_beta_"), m))
      r1 <- na3("OracleG")
      r2 <- score_test_ObsG_single(G_all, step1_pars = ObsG_typeIerror_step1(mydf))
      r3 <- score_test_SynSurrG_single(G_all, step1_pars = SynSurrG_typeIerror_step1(mydf))
      r4 <- na3("Oracle")
      r5 <- na3("Obs")
      r6 <- score_test_SynSurr_single(G_all, step1_pars = SynSurr_typeIerror_step1(mydf, independent_indices), independent_indices)
      res[ss, ] <- c(r1$negative_log10_pval_OracleG, r1$hat_beta_OracleG, r1$var_hat_beta_OracleG,
                     r2$negative_log10_pval_ObsG, r2$hat_beta_ObsG, r2$var_hat_beta_ObsG,
                     r3$negative_log10_pval_SynSurrG, r3$hat_beta_SynSurrG, r3$var_hat_beta_SynSurrG,
                     r4$negative_log10_pval_Oracle, r4$hat_beta_Oracle, r4$var_hat_beta_Oracle,
                     r5$negative_log10_pval_Obs, r5$hat_beta_Obs, r5$var_hat_beta_Obs,
                     r6$negative_log10_pval_SynSurr, r6$hat_beta_SynSurr, r6$var_hat_beta_SynSurr)
    }
    attr(res, "snps") <- rs_array
    attr(res, "methods") <- c("ObsG", "SynSurrG", "SynSurr")
    attr(res, "settings") <- list(miss = miss, rho = rho, pve_g = pve_g, n_obs = n_obs, n_miss = n_miss,
                                  array = array_id, arraynum = arraynum, seed = seed,
                                  seed_fix_constant = seed_g, snp_sample_seed = 241150L,
                                  cov_TS = sim_cov_TS(rho), genotypes = paste0(geno_prefix, " (autosomes)"))
    saveRDS(res, out_file(miss, pve_g))
    tick("miss ", miss, " pve_g ", pve_g, ": saved (seed ", seed, ")")
  }
  rm(pre_vars_step2, GRM, Genotype_matrix); invisible(gc())
}

## SUPERSEDED (2026-10-14) by true_geno_typeIerror_calgeno_rev1.R (all four rho, 17 replicates); not run in production.
## Rerun of the rho = 0.75 settings of the true-genotype type I error simulation
## (Supp Table "typeIcontrol_ture"; simluation_typeIerror_control_true_geno_v2.R),
## whose genetic covariance broke Cauchy-Schwarz (Reviewer 2, comment 5).
## Same data generation (UKB sparse GRM, the first n = n_obs + n_miss eligible
## White-British eids in the shuffled GRM order, DGP_GRMr_nullmodel_step1-3),
## same methods and output columns as the submitted script; changes:
##   - rho = 0.75 only, with sim_cov_TS() (genetic correlation 0.99, total rho 0.75);
##   - autosomal genotyped SNPs only (chromosomes 1-22 of
##     genotype_data/ukb_cal_allchr_v2_whiteb_QC.bim; 575,236 of 591,558);
##   - every random draw is seeded: one seed per (array, missing rate), set before
##     the null phenotype is drawn (seed saved in attr(, "settings"));
##   - Syn-PALM is scored with score_test_SynSurrG_fast (funcs_fastscan.R; identical
##     to score_test_SynSurrG_multiply up to 1e-13); the other five tests use the
##     *_multiply functions of the submitted script;
##   - one output file per (missing rate, array), skipped if it already exists, so
##     a preempted / requeued task resumes where it stopped.
## As in the submitted script, every array draws its own null phenotype, fits the
## six null models once and scores its block of SNPs in chunks of 200.
##   Rscript true_geno_typeIerror_rho075_rev1.R <array> [n_arrays=1000] [miss list, e.g. 0,0.5] [max_snps] [out_dir]
## Output: simulation_code/results_rev1/true_geno_typeIerror_rho075/
##   missing_<m>_rho_0.75_chunkarray_<a>_negative_log10_pvals_type1error_rev1.rds
##   (18 columns as the submitted temp files; rsids in attr(, "snps")).
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R") })
options(synsurrg.verbose = FALSE)
library(withr)

args     <- commandArgs(trailingOnly = TRUE)
array_id <- as.integer(args[1])
arraynum <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 1000L
miss_all <- c(0, 0.25, 0.5, 0.75, 0.9)
miss_vec <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(strsplit(args[3], ",")[[1]]) else miss_all
max_snps <- if (length(args) >= 4 && nzchar(args[4])) as.integer(args[4]) else NA_integer_
OUT      <- if (length(args) >= 5 && nzchar(args[5])) args[5] else "simulation_code/results_rev1/true_geno_typeIerror_rho075/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rho <- 0.75
n_obs <- 3e4
tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5
pve_x <- 0
chunk_size <- 200

out_file <- function(miss) sprintf("%smissing_%s_rho_%s_chunkarray_%d_negative_log10_pvals_type1error_rev1.rds",
                                   OUT, miss, rho, array_id)
todo <- miss_vec[!file.exists(vapply(miss_vec, out_file, ""))]
if (!length(todo)) { cat("all outputs exist for array", array_id, "\n"); quit(save = "no") }

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
rs_auto <- colnames(G0)[bim$V1 %in% as.character(1:22)]          # autosomes only
rs_index_list <- split(rs_auto, cut(seq_along(rs_auto), breaks = arraynum, labels = FALSE))
rs_array <- rs_index_list[[array_id]]
if (!is.na(max_snps)) rs_array <- head(rs_array, max_snps)
tick("array ", array_id, "/", arraynum, ": ", length(rs_array), " autosomal SNPs, miss ", paste(todo, collapse = ","))

for (miss in todo) {
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

  seed <- 20261014L + 100L * array_id + mi
  set.seed(seed)
  pre_vars_step1 <- DGP_GRMr_nullmodel_step1(n_obs, miss, tauT, tauS, sigmaT, sigmaS, GRM)
  pre_vars_step2 <- DGP_GRMr_nullmodel_step2(n_obs, miss, tauT, tauS, sigmaT, sigmaS, rho, GRM, pre_vars_step1)
  rm(pre_vars_step1); invisible(gc())
  mydf <- DGP_GRMr_nullmodel_step3(n_obs, miss, tauT, sigmaT, pve_x, pre_vars_step2, GRM)
  rm(pre_vars_step2); invisible(gc())

  SynSurrG_pars <- add_fastscan_cache(SynSurrG_typeIerror_step1(mydf))
  OracleG_pars  <- OracleG_typeIerror_step1(mydf)
  ObsG_pars     <- ObsG_typeIerror_step1(mydf)
  SynSurr_pars  <- SynSurr_typeIerror_step1(mydf, independent_indices)
  Oracle_pars   <- Oracle_typeIerror_step1(mydf, independent_indices)
  Obs_pars      <- Obs_typeIerror_step1(mydf, independent_indices)
  tick("miss ", miss, ": null models fitted (n = ", n_all, ", seed ", seed, ")")

  matched_indices <- match(subeids[1:n_all], geno_samp.id)
  if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")

  num_chunks <- ceiling(length(rs_array) / chunk_size)
  res <- vector("list", num_chunks)
  for (chunk_index in seq_len(num_chunks)) {
    chunk_idx <- rs_array[((chunk_index - 1) * chunk_size + 1):min(chunk_index * chunk_size, length(rs_array))]
    Genotype_matrix <- as.matrix(G0[matched_indices, chunk_idx, drop = FALSE])
    Genotype_matrix <- replace(Genotype_matrix, is.na(Genotype_matrix), 0)
    Genotype_matrix <- fix_constant_columns(Genotype_matrix, intersect(1:n_obs, independent_indices))

    r_sg <- score_test_SynSurrG_fast(Genotype_matrix, SynSurrG_pars)
    r_og <- score_test_OracleG_multiply(g_matrix = Genotype_matrix, step1_pars = OracleG_pars)
    r_bg <- score_test_ObsG_multiply(g_matrix = Genotype_matrix, step1_pars = ObsG_pars)
    r_s  <- score_test_SynSurr_multiply(g_matrix = Genotype_matrix, step1_pars = SynSurr_pars, independent_indices)
    r_o  <- score_test_Oracle_multiply(g_matrix = Genotype_matrix, step1_pars = Oracle_pars, independent_indices)
    r_b  <- score_test_Obs_multiply(g_matrix = Genotype_matrix, step1_pars = Obs_pars, independent_indices)

    res[[chunk_index]] <- cbind(OracleG_negative_log10_pvalues  = r_og$negative_log10_pval_OracleG,
                                OracleG_hat_beta                = r_og$hat_beta_OracleG,
                                OracleG_var_hat_beta            = r_og$var_hat_beta_OracleG,
                                SynSurrG_negative_log10_pvalues = r_sg$negative_log10_pval_SynSurrG,
                                SynSurrG_hat_beta               = r_sg$hat_beta_SynSurrG,
                                SynSurrG_var_hat_beta           = r_sg$var_hat_beta_SynSurrG,
                                ObsG_negative_log10_pvalues     = r_bg$negative_log10_pval_ObsG,
                                ObsG_hat_beta                   = r_bg$hat_beta_ObsG,
                                ObsG_var_hat_beta               = r_bg$var_hat_beta_ObsG,
                                Oracle_negative_log10_pvalues   = r_o$negative_log10_pval_Oracle,
                                Oracle_hat_beta                 = r_o$hat_beta_Oracle,
                                Oracle_var_hat_beta             = r_o$var_hat_beta_Oracle,
                                SynSurr_negative_log10_pvalues  = r_s$negative_log10_pval_SynSurr,
                                SynSurr_hat_beta                = r_s$hat_beta_SynSurr,
                                SynSurr_var_hat_beta            = r_s$var_hat_beta_SynSurr,
                                Obs_negative_log10_pvalues      = r_b$negative_log10_pval_Obs,
                                Obs_hat_beta                    = r_b$hat_beta_Obs,
                                Obs_var_hat_beta                = r_b$var_hat_beta_Obs)
    tick("miss ", miss, ": chunk ", chunk_index, "/", num_chunks)
  }
  results <- do.call(rbind, res)
  stopifnot(nrow(results) == length(rs_array))
  attr(results, "snps") <- rs_array
  attr(results, "settings") <- list(miss = miss, rho = rho, n_obs = n_obs, n_miss = n_miss, array = array_id,
                                    arraynum = arraynum, seed = seed, cov_TS = sim_cov_TS(rho),
                                    genotypes = paste0(geno_prefix, " (autosomes)"),
                                    vc = SynSurrG_pars$vc)
  saveRDS(results, out_file(miss))
  tick("saved ", out_file(miss))
  rm(mydf, SynSurrG_pars, OracleG_pars, ObsG_pars, SynSurr_pars, Oracle_pars, Obs_pars, res, results, GRM)
  invisible(gc())
}

## True-genotype type I error simulation of the revision (Supp Table
## "typeIcontrol_ture"), all four rho on the calibrated (genotyped) autosomal SNPs.
## Replaces the submitted imputed-variant table, whose generating script and
## hard-call genotype files no longer exist. Design as
## simluation_typeIerror_control_true_geno_v2.R / true_geno_typeIerror_rho075_rev1.R:
##   - UKB sparse GRM (blocks < 500, shuffled with seed 255155), the first
##     n = n_obs + n_miss eligible White-British eids, independent_indices seed 25324;
##   - DGP_GRMr_nullmodel_step1-3 of funcs_v2.R with sim_cov_TS() for every rho
##     (rho = 0.75: genetic correlation 0.99, total rho 0.75);
##   - genotypes: ukb_cal_allchr_v2_whiteb_QC, chromosomes 1-22 (575,236 SNPs),
##     NA -> 0, fix_constant_columns on the observed independent individuals;
##   - all autosomal SNPs are scored under R null phenotype replicates; the SNPs are
##     split into n_arrays blocks and every (replicate, array, missing rate, rho)
##     draws its own null phenotype and fits its own null models;
##   - only Syn-PALM (score_test_SynSurrG_fast after add_fastscan_cache; identical to
##     score_test_SynSurrG_multiply up to 1e-13) and Obs-GWAS (score_test_ObsG_multiply)
##     are computed; Oracle / SynSurr / Obs are not needed for the table.
## Seeds (never collide for replicate <= 99, array <= 999):
##   phenotype  seed = 30261015 + 100000 * replicate + 100 * array + 10 * miss_index + rho_index
##   fix_constant_columns seed = same with rho_index = 0
## A chunk of 200 SNPs is read once per missing rate and scored against the four rho.
## One file per (replicate, array, missing rate, rho); a (replicate, array, missing
## rate) whose four files exist is skipped, so requeued / resubmitted tasks resume.
##   Rscript true_geno_typeIerror_calgeno_rev1.R <replicate> <array> [n_arrays=100] [miss list] [max_snps] [out_dir]
## Output: simulation_code/results_rev1/true_geno_typeIerror_calgeno/
##   rep_<r>_missing_<m>_rho_<rho>_chunkarray_<a>_negative_log10_pvals_type1error_rev1.rds
##   columns SynSurrG_{negative_log10_pvalues,hat_beta,var_hat_beta}, ObsG_{...}
##   (names as in the submitted files); rsids in attr(, "snps"), seeds in attr(, "settings").
setwd("<PROJECT_DIR>")
suppressMessages({ source("code/funcs_v2.R"); source("code/fastscan/funcs_fastscan.R") })
options(synsurrg.verbose = FALSE)
library(withr)

args      <- commandArgs(trailingOnly = TRUE)
replicate <- as.integer(args[1])
array_id  <- as.integer(args[2])
arraynum  <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 100L
miss_all  <- c(0, 0.25, 0.5, 0.75, 0.9)
miss_vec  <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(strsplit(args[4], ",")[[1]]) else miss_all
max_snps  <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else NA_integer_
OUT       <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "simulation_code/results_rev1/true_geno_typeIerror_calgeno/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(replicate >= 1, replicate <= 99, array_id >= 1, array_id <= arraynum, arraynum <= 999)

rho_vec <- c(0, 0.25, 0.5, 0.75)
n_obs <- 3e4
tauT <- 0.7; tauS <- 0.4; sigmaT <- 0.7; sigmaS <- 0.5
pve_x <- 0
chunk_size <- 200
seed_of <- function(mi, ri) 30261015L + 100000L * replicate + 100L * array_id + 10L * mi + ri

out_file <- function(miss, rho) sprintf("%srep_%d_missing_%s_rho_%s_chunkarray_%d_negative_log10_pvals_type1error_rev1.rds",
                                        OUT, replicate, miss, rho, array_id)
done <- vapply(miss_vec, function(m) all(file.exists(vapply(rho_vec, function(r) out_file(m, r), ""))), TRUE)
todo <- miss_vec[!done]
if (!length(todo)) { cat("all outputs exist for replicate", replicate, "array", array_id, "\n"); quit(save = "no") }

t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%7.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

### GRM, sample and genotypes: as in the submitted v2 script ################
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
rs_array <- split(rs_auto, cut(seq_along(rs_auto), breaks = arraynum, labels = FALSE))[[array_id]]
if (!is.na(max_snps)) rs_array <- head(rs_array, max_snps)
tick("replicate ", replicate, ", array ", array_id, "/", arraynum, ": ", length(rs_array),
     " autosomal SNPs, miss ", paste(todo, collapse = ","))

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
  pre_vars_step1 <- DGP_GRMr_nullmodel_step1(n_obs, miss, tauT, tauS, sigmaT, sigmaS, GRM)

  ## one null phenotype per rho, null models fitted once
  fits <- lapply(seq_along(rho_vec), function(ri) {
    pre_vars_step2 <- DGP_GRMr_nullmodel_step2(n_obs, miss, tauT, tauS, sigmaT, sigmaS, rho_vec[ri], GRM, pre_vars_step1)
    set.seed(seed_of(mi, ri))
    mydf <- DGP_GRMr_nullmodel_step3(n_obs, miss, tauT, sigmaT, pve_x, pre_vars_step2, GRM)
    list(syn = add_fastscan_cache(SynSurrG_typeIerror_step1(mydf)), obs = ObsG_typeIerror_step1(mydf))
  })
  rm(pre_vars_step1); invisible(gc())
  tick("miss ", miss, ": null models fitted (n = ", n_all, ")")

  matched_indices <- match(subeids[1:n_all], geno_samp.id)
  if (any(is.na(matched_indices))) stop("Some eids in test_df are not found in G.")
  set.seed(seed_of(mi, 0L))                                     # fix_constant_columns draws
  num_chunks <- ceiling(length(rs_array) / chunk_size)
  res <- vector("list", length(rho_vec))
  for (chunk_index in seq_len(num_chunks)) {
    chunk_idx <- rs_array[((chunk_index - 1) * chunk_size + 1):min(chunk_index * chunk_size, length(rs_array))]
    Genotype_matrix <- as.matrix(G0[matched_indices, chunk_idx, drop = FALSE])
    Genotype_matrix <- replace(Genotype_matrix, is.na(Genotype_matrix), 0)
    Genotype_matrix <- fix_constant_columns(Genotype_matrix, intersect(1:n_obs, independent_indices))
    for (ri in seq_along(rho_vec)) {
      r_sg <- score_test_SynSurrG_fast(Genotype_matrix, fits[[ri]]$syn)
      r_bg <- score_test_ObsG_multiply(g_matrix = Genotype_matrix, step1_pars = fits[[ri]]$obs)
      res[[ri]][[chunk_index]] <- cbind(SynSurrG_negative_log10_pvalues = r_sg$negative_log10_pval_SynSurrG,
                                        SynSurrG_hat_beta               = r_sg$hat_beta_SynSurrG,
                                        SynSurrG_var_hat_beta           = r_sg$var_hat_beta_SynSurrG,
                                        ObsG_negative_log10_pvalues     = r_bg$negative_log10_pval_ObsG,
                                        ObsG_hat_beta                   = r_bg$hat_beta_ObsG,
                                        ObsG_var_hat_beta               = r_bg$var_hat_beta_ObsG)
    }
    if (chunk_index %% 5 == 0 || chunk_index == num_chunks) tick("miss ", miss, ": chunk ", chunk_index, "/", num_chunks)
  }
  for (ri in seq_along(rho_vec)) {
    results <- do.call(rbind, res[[ri]])
    stopifnot(nrow(results) == length(rs_array))
    attr(results, "snps") <- rs_array
    attr(results, "settings") <- list(replicate = replicate, array = array_id, arraynum = arraynum,
                                      miss = miss, rho = rho_vec[ri], n_obs = n_obs, n_miss = n_miss,
                                      seed = seed_of(mi, ri), seed_fix_constant = seed_of(mi, 0L),
                                      cov_TS = sim_cov_TS(rho_vec[ri]),
                                      genotypes = paste0(geno_prefix, " (autosomes)"),
                                      vc = fits[[ri]]$syn$vc)
    saveRDS(results, out_file(miss, rho_vec[ri]))
  }
  tick("miss ", miss, ": saved 4 rho files")
  rm(fits, res, GRM); invisible(gc())
}

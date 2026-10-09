## ---------------------------------------------------------------------------
## SynPALM on a real biobank (e.g. All of Us): template
##
## Two stages, so the genome-wide scan can run as many parallel jobs:
##
##   Rscript biobank_template.R fit  <protein>          # once per protein
##   Rscript biobank_template.R scan <protein> <chr>    # once per chromosome
##
## Stage "fit" builds the folds, the 5-fold cross-fitted random-forest
## synthetic phenotype and the null model, and saves them; it also writes
## <protein>_rf_accuracy.tsv, the correlation between synthetic and measured
## protein per fold and overall. Stage "scan" reads that file
## and tests one chromosome. Every path below is a placeholder: point them at
## your own files. Sample IDs must be the same strings in every input.
## ---------------------------------------------------------------------------

library(SynPALM)
library(Matrix)
library(data.table)

## --- paths (placeholders) --------------------------------------------------
PROTEIN_FILE    <- "data/proteins.tsv"         # id + one column per protein; may list measured people only
SURROGATE_FILE  <- "data/surrogates.tsv"       # id + candidate surrogates (labs, vitals, ...)
DIAGNOSIS_FILE  <- "data/diagnoses.tsv"        # id + 0/1 indicators (optional; set NULL)
COVARIATE_FILE  <- "data/covariates.tsv"       # id + age, sex, PCs, batch, ...
GRM_FILE        <- "data/sparse_grm.rds"       # sparse GRM, dimnames = sample IDs
GENOTYPE_PREFIX <- "genotypes/chr%s.bed"       # PLINK .bed per chromosome (+ .bim/.fam)
OUT_DIR         <- "results"

## --- analysis settings -----------------------------------------------------
RF_FIXED      <- c("age", "sex")   # surrogate columns always given to the forest
N_TOP         <- 100               # surrogates kept by per-fold screening
METHODS       <- c("SynSurrG", "ObsG")
VC_CONSTRAINT <- "none"            # see ?SynPALM-package, "Numerical stability"
THREADS       <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

options(synsurrg.ncores = THREADS)

args  <- commandArgs(trailingOnly = TRUE)
stage <- args[1]
this_protein <- args[2]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
fit_file <- file.path(OUT_DIR, paste0(this_protein, "_synpalm_fit.rds"))

read_by_id <- function(path) {
  x <- fread(path, data.table = FALSE)
  rownames(x) <- as.character(x[[1]])
  x[-1]
}

if (stage == "fit") {
  prot <- fread(PROTEIN_FILE, select = c(1, which(names(fread(PROTEIN_FILE, nrows = 0)) == this_protein)),
                data.table = FALSE)
  protein     <- setNames(prot[[2]], as.character(prot[[1]]))
  rf_features <- read_by_id(SURROGATE_FILE)
  rf_binary   <- if (!is.null(DIAGNOSIS_FILE)) read_by_id(DIAGNOSIS_FILE)
  covariates  <- read_by_id(COVARIATE_FILE)
  grm         <- readRDS(GRM_FILE)

  ## the analysis set: everyone in the GRM, surrogates and covariates (in GRM
  ## order), measured or not. It is NOT restricted to the protein table, which
  ## usually lists only the measured people: the unmeasured majority is what
  ## the synthetic phenotype adds. Absent from the protein table = unmeasured.
  keep <- Reduce(intersect, list(rownames(rf_features), rownames(covariates),
                                 if (!is.null(rf_binary)) rownames(rf_binary)))
  ids  <- rownames(grm)[rownames(grm) %in% keep]
  grm  <- grm[ids, ids]
  protein <- setNames(protein[match(ids, names(protein))], ids)
  cat(this_protein, ": ", length(ids), " individuals, ",
      sum(!is.na(protein)), " with a measured protein\n", sep = "")

  folds <- synpalm_folds(grm, K = 5)
  pred  <- synpalm_predict(protein, rf_features, folds, rf_binary = rf_binary,
                           rf_fixed = RF_FIXED, n_top = N_TOP, num.threads = THREADS)
  null  <- synpalm_null(pred$protein, pred$protein_hat, covariates, grm, folds,
                        methods = METHODS, vc_constraint = VC_CONSTRAINT)

  ## accuracy of the synthetic phenotype (5-fold cross-fit): one row per fold + "all"
  acc <- cbind(protein = this_protein, pred$accuracy)
  print(acc, digits = 4, row.names = FALSE)
  fwrite(acc, file.path(OUT_DIR, paste0(this_protein, "_rf_accuracy.tsv")), sep = "\t")

  vc <- null$pars$SynSurrG$vc
  cat("variance components valid without projection:",
      !(isTRUE(vc$would_project_G) || isTRUE(vc$would_project_E)), "\n")
  saveRDS(list(null = null, prediction = pred[c("ids", "protein_hat", "rho_oof", "rho_by_fold",
                                                 "accuracy", "n_obs", "fold_features")]),
          fit_file)
  cat("written:", fit_file, "\n")

} else if (stage == "scan") {
  chr  <- args[3]
  fit  <- readRDS(fit_file)
  res  <- synpalm_scan(fit$null, sprintf(GENOTYPE_PREFIX, chr))
  res$protein <- this_protein
  res$rho_oof <- fit$prediction$rho_oof   # synthetic vs measured protein, cross-fitted
  out  <- file.path(OUT_DIR, sprintf("%s_chr%s_synpalm.tsv.gz", this_protein, chr))
  fwrite(res, out, sep = "\t")
  cat("written: ", out, " (", nrow(res), " variants)\n", sep = "")

} else {
  stop("usage: Rscript biobank_template.R fit <protein> | scan <protein> <chr>")
}

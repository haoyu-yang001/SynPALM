## ---------------------------------------------------------------------------
## SynPALM quick-start demo
##
## Self-contained: simulates a cohort with family structure, a protein measured
## in only 10% of it, many candidate surrogate variables for the random forest,
## covariates and genotypes, so no external data are needed. Then runs the
## whole analysis with one call:
##
##   synpalm_gwas()  =  folds -> cross-fitted RF -> null model -> score tests
##
## Run with:
##   source(system.file("examples", "quickstart.R", package = "SynPALM"))
##
## The seed is fixed so the printed output matches the "Expected output"
## section of README.md.
## ---------------------------------------------------------------------------

library(SynPALM)
library(Matrix)

## --- settings --------------------------------------------------------------
N_FAM     <- 5000    # number of families
FAM_SIZE  <- 4       # relatives per family  -> cohort = N_FAM * FAM_SIZE
MISS_RATE <- 0.90    # fraction of the cohort with NO measured protein
N_SNP     <- 500     # variants tested; the first one is causal
BETA_G    <- 0.08    # causal effect size, standardised genotype scale
TAU       <- 0.40    # polygenic variance component
SIGMA     <- 0.60    # residual variance component
N_SURR    <- 200     # candidate surrogate variables offered to the forest
N_INFORM  <- 10      # of which this many carry information on the protein
N_DISEASE <- 30      # binary (disease-type) candidate surrogates

set.seed(1)
demo_start <- Sys.time()
n_all <- N_FAM * FAM_SIZE
ids   <- sprintf("ID%06d", seq_len(n_all))

## ---------------------------------------------------------------------------
## 1. Sparse GRM with family structure (within-family kinship 0.5)
## ---------------------------------------------------------------------------
fam_block <- matrix(0.5, FAM_SIZE, FAM_SIZE)
diag(fam_block) <- 1
GRM <- Matrix::bdiag(replicate(N_FAM, fam_block, simplify = FALSE))
GRM <- methods::as(methods::as(GRM, "CsparseMatrix"), "generalMatrix")
dimnames(GRM) <- list(ids, ids)

## ---------------------------------------------------------------------------
## 2. Covariates for the mixed model
## ---------------------------------------------------------------------------
covariates <- data.frame(
  age = as.numeric(scale(round(rnorm(n_all, 57, 8)))),
  sex = factor(rbinom(n_all, 1, 0.5), labels = c("F", "M")),
  matrix(rnorm(n_all * 5), n_all, 5, dimnames = list(NULL, paste0("PC", 1:5))),
  row.names = ids
)

## ---------------------------------------------------------------------------
## 3. Genotypes (allele counts 0/1/2), the first variant causal
## ---------------------------------------------------------------------------
maf  <- runif(N_SNP, 0.2, 0.4)
Gmat <- vapply(seq_len(N_SNP), function(j) rbinom(n_all, 2, maf[j]), numeric(n_all))
dimnames(Gmat) <- list(ids, paste0("rs_demo_", seq_len(N_SNP)))

## ---------------------------------------------------------------------------
## 4. Protein under the mixed model:  y = X b + g + u + e,  u ~ N(0, TAU * GRM)
## ---------------------------------------------------------------------------
L <- Matrix::t(Matrix::chol(GRM))
u <- sqrt(TAU) * as.vector(L %*% rnorm(n_all))
y <- 0.3 * covariates$age + 0.2 * (covariates$sex == "M") +
  BETA_G * as.numeric(scale(Gmat[, 1])) + u + rnorm(n_all, 0, sqrt(SIGMA))

## ---------------------------------------------------------------------------
## 5. Candidate surrogates for the random forest
##
## Stand-ins for clinical measurements and diagnoses. Only N_INFORM of the
## N_SURR continuous ones, and a few of the binary ones, relate to the protein;
## the forest's per-fold screening has to find them.
## ---------------------------------------------------------------------------
w  <- c(seq(0.9, 0.3, length.out = N_INFORM), rep(0, N_SURR - N_INFORM))
rf_features <- as.data.frame(sapply(w, function(a) a * y + rnorm(n_all)))
names(rf_features) <- paste0("lab", seq_len(N_SURR))
rf_features$age <- covariates$age
rf_features$sex <- as.numeric(covariates$sex == "M")
rownames(rf_features) <- ids

rf_binary <- as.data.frame(sapply(seq_len(N_DISEASE), function(k)
  rbinom(n_all, 1, plogis(-1.5 + if (k <= 3) 0.8 * y else 0))))
names(rf_binary) <- paste0("dx", seq_len(N_DISEASE))
rownames(rf_binary) <- ids

## ---------------------------------------------------------------------------
## 6. The protein is measured in only (1 - MISS_RATE) of the cohort
## ---------------------------------------------------------------------------
protein <- setNames(rep(NA_real_, n_all), ids)
obs <- sort(sample.int(n_all, floor((1 - MISS_RATE) * n_all)))
protein[obs] <- y[obs]
cat("Cohort: ", n_all, " individuals | ", length(obs), " with a measured protein | ",
    N_SURR, " continuous + ", N_DISEASE, " binary candidate surrogates\n\n", sep = "")

## ---------------------------------------------------------------------------
## 7. The analysis
## ---------------------------------------------------------------------------
options(synsurrg.verbose = FALSE)
t_run <- system.time(
  fit <- synpalm_gwas(protein     = protein,
                      covariates  = covariates,
                      rf_features = rf_features,
                      rf_binary   = rf_binary,
                      rf_fixed    = c("age", "sex"),
                      grm         = GRM,
                      genotype    = Gmat,
                      verbose     = FALSE)
)

## ---------------------------------------------------------------------------
## 8. Results
## ---------------------------------------------------------------------------
pr <- fit$prediction
cat("--- Synthetic phenotype (5-fold cross-fit) ---\n")
cat("  correlation with the measured protein, per fold and overall:\n")
print(fit$accuracy, digits = 3, row.names = FALSE)
sel <- table(unlist(lapply(pr$fold_features, `[[`, "cov_pred")))
cat("  informative surrogates selected in all 5 folds   : ",
    sum(sel[intersect(names(sel), paste0("lab", seq_len(N_INFORM)))] == 5), " of ", N_INFORM, "\n", sep = "")

res <- fit$results
cat("\n--- First rows of fit$results ---\n")
print(head(res, 3), digits = 3)

cat("\n--- Causal variant (true standardised effect = ", BETA_G, ") ---\n", sep = "")
print(data.frame(method = c("SynPALM (SynSurrG)", "observed only (ObsG)"),
                 log10p = round(c(res$SynSurrG_log10p[1], res$ObsG_log10p[1]), 3),
                 beta   = round(c(res$SynSurrG_beta[1], res$ObsG_beta[1]), 4),
                 se     = round(c(res$SynSurrG_se[1], res$ObsG_se[1]), 4)),
      row.names = FALSE)

cat("\n--- Calibration on ", N_SNP - 1, " null variants ---\n", sep = "")
for (m in c("SynSurrG", "ObsG")) {
  p <- res[[paste0(m, "_p")]][-1]
  lambda <- median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  cat(sprintf("  %-8s lambda = %.3f   P<0.05: %.3f   P<0.01: %.3f\n",
              m, lambda, mean(p < 0.05), mean(p < 0.01)))
}

vc <- fit$null$pars$SynSurrG$vc
cat("\n--- Variance components (SynSurrG) ---\n")
print(round(vc$used, 4))
cat("  valid without projection: ", !(vc$would_project_G || vc$would_project_E), "\n", sep = "")

cat("\n--- Timing ---\n")
cat("synpalm_gwas() elapsed (s): ", round(t_run[["elapsed"]], 1), "\n", sep = "")
cat("Whole demo elapsed (s)    : ",
    round(as.numeric(difftime(Sys.time(), demo_start, units = "secs")), 1), "\n", sep = "")
cat(R.version.string, "|", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")

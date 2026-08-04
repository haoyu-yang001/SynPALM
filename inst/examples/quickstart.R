## ---------------------------------------------------------------------------
## SynPALM quick-start demo
##
## Self-contained: simulates a cohort with family structure and 90% missingness
## in the target phenotype, so no external data are needed. It follows the same
## call sequence used for the UK Biobank analysis in the manuscript:
##
##   build mydf  ->  SynSurrG_ablation_estimate()  ->  score_test_SynSurrG_multiply()
##
## Run with:
##   source(system.file("examples", "quickstart.R", package = "SynPALM"))
##
## The seed is fixed so the printed output matches the "Expected output"
## section of README.md. Do not change it.
##
## Scale: with the settings below the demo needs roughly 2-4 GB of RAM. To run
## it on a smaller machine, reduce N_FAM (cohort size) and N_SNP; both scale
## the cost roughly linearly.
## ---------------------------------------------------------------------------

library(SynPALM)
library(Matrix)

## --- settings --------------------------------------------------------------
N_FAM     <- 5000    # number of families
FAM_SIZE  <- 4       # relatives per family  -> cohort = N_FAM * FAM_SIZE
MISS_RATE <- 0.90    # fraction of the cohort with NO measured target phenotype
N_SNP     <- 500     # variants tested; the first one is causal
BETA_G    <- 0.08    # causal effect size, standardised scale
R_TS      <- 0.60    # correlation between target and synthetic phenotype
TAU       <- 0.40    # polygenic variance component
SIGMA     <- 0.60    # residual variance component

set.seed(1)
demo_start <- Sys.time()

n_all    <- N_FAM * FAM_SIZE
obs_frac <- 1 - MISS_RATE

## ---------------------------------------------------------------------------
## 1. A sparse GRM with family structure
##
## SynPALM targets cohorts with cryptic relatedness, so the demo uses genuine
## family blocks rather than an identity matrix. Within-family kinship 0.5.
## ---------------------------------------------------------------------------
fam_block <- matrix(0.5, FAM_SIZE, FAM_SIZE)
diag(fam_block) <- 1
GRM <- Matrix::bdiag(replicate(N_FAM, fam_block, simplify = FALSE))
GRM <- methods::as(methods::as(GRM, "CsparseMatrix"), "generalMatrix")

cat("GRM: ", nrow(GRM), " x ", ncol(GRM),
    "  | class: ", class(GRM)[1],
    "  | nonzero fraction: ",
    signif(Matrix::nnzero(GRM) / prod(dim(GRM)), 3), "\n", sep = "")

## ---------------------------------------------------------------------------
## 2. Covariates
## ---------------------------------------------------------------------------
X_all <- data.frame(
  age = scale(round(rnorm(n_all, 57, 8)))[, 1],
  sex = rbinom(n_all, 1, 0.5),
  matrix(rnorm(n_all * 5), n_all, 5)
)
colnames(X_all) <- c("age", "sex", paste0("geneticPC", 1:5))

## ---------------------------------------------------------------------------
## 3. Genotypes -- N_SNP variants, the first of which is causal
## ---------------------------------------------------------------------------
maf  <- runif(N_SNP, 0.2, 0.4)
Gmat <- vapply(seq_len(N_SNP), function(j) rbinom(n_all, 2, maf[j]), numeric(n_all))
colnames(Gmat) <- paste0("rs_demo_", seq_len(N_SNP))

## ---------------------------------------------------------------------------
## 4. Target phenotype under the mixed model
##      Y = X beta + g_causal + u + e ,   u ~ N(0, TAU * GRM)
## ---------------------------------------------------------------------------
L    <- Matrix::t(Matrix::chol(GRM))
u    <- sqrt(TAU) * as.vector(L %*% rnorm(n_all))
beta <- c(0.3, 0.2, rep(0.05, 5))

Y_true <- as.vector(as.matrix(X_all) %*% beta) +
  BETA_G * scale(Gmat[, 1])[, 1] +
  u +
  rnorm(n_all, 0, sqrt(SIGMA))

## ---------------------------------------------------------------------------
## 5. Synthetic / surrogate phenotype
##
## Stands in for the machine-learning prediction (a random forest on clinical
## covariates in the manuscript). SynPALM gains power as its correlation with
## the target rises, and stays valid when that correlation is poor.
## ---------------------------------------------------------------------------
int_transform <- function(v) {
  r <- rank(v, na.last = "keep")
  n <- sum(!is.na(v))
  qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
}

S <- int_transform(R_TS * scale(Y_true)[, 1] + sqrt(1 - R_TS^2) * rnorm(n_all))

## ---------------------------------------------------------------------------
## 6. Missingness
## ---------------------------------------------------------------------------
obs_idx <- sort(sample.int(n_all, floor(obs_frac * n_all)))

Y_obs <- rep(NA_real_, n_all)
Y_obs[obs_idx] <- int_transform(Y_true[obs_idx])

cat("Cohort: ", n_all, " individuals | ", length(obs_idx),
    " with an observed target phenotype (", round(100 * obs_frac), "%), ",
    "missingness ", round(100 * MISS_RATE), "%\n", sep = "")
cat("Observed-vs-synthetic correlation: ",
    round(cor(Y_obs[obs_idx], S[obs_idx]), 3), "\n\n", sep = "")

## ---------------------------------------------------------------------------
## 7. Assemble mydf -- the single input object SynPALM entry points take
## ---------------------------------------------------------------------------
mydf <- list(
  X_all = X_all,   # covariates for all individuals, no intercept column
  S     = S,       # synthetic phenotype, complete
  Y_obs = Y_obs,   # target phenotype, NA where unmeasured
  GRM   = GRM      # sparse GRM over all individuals
)

cat("--- structure of mydf ---\n")
str(mydf, max.level = 1, give.attr = FALSE)

## ---------------------------------------------------------------------------
## 8. Relatedness blocks
##
## Also a check on block detection: with N_FAM disjoint families the algorithm
## must recover exactly N_FAM blocks.
## ---------------------------------------------------------------------------
blocks <- find_blocks_vectorized(GRM)
cat("\nRelatedness blocks found: ", length(blocks),
    " (expected ", N_FAM, ")\n", sep = "")
stopifnot(length(blocks) == N_FAM)

independent_indices <- vapply(blocks,
                              function(b) b[sample.int(length(b), 1)],
                              integer(1))
obs_protein_index <- which(!is.na(mydf$Y_obs))

Gmat <- fix_constant_columns(Gmat,
                             intersect(obs_protein_index, independent_indices))

## ---------------------------------------------------------------------------
## 9. Step 1 -- null model / variance components
## ---------------------------------------------------------------------------
cat("\nFitting the null model (step 1)...\n")
t_step1 <- system.time(
  SynSurrG_step1 <- SynSurrG_ablation_estimate(mydf)
)
cat("Step 1 elapsed (s): ", round(t_step1[["elapsed"]], 1), "\n", sep = "")

## ---------------------------------------------------------------------------
## 10. Step 2 -- score test across all variants
## ---------------------------------------------------------------------------
cat("\nRunning the SynPALM score test (step 2)...\n")
t_test <- system.time(
  res <- score_test_SynSurrG_multiply(g_matrix = Gmat, step1_pars = SynSurrG_step1)
)
cat("Score test elapsed (s): ", round(t_test[["elapsed"]], 1),
    "  (", N_SNP, " variants, ",
    round(1000 * t_test[["elapsed"]] / N_SNP, 1), " ms per variant)\n", sep = "")

## ---------------------------------------------------------------------------
## 11. Results
## ---------------------------------------------------------------------------
cat("\n--- components returned by score_test_SynSurrG_multiply ---\n")
print(names(res))

nl10  <- as.numeric(res$negative_log10_pval_SynSurrG)
tstat <- as.numeric(res$T_score_SynSurrG)
bhat  <- as.numeric(res$hat_beta_SynSurrG)
vhat  <- as.numeric(res$var_hat_beta_SynSurrG)

cat("\n--- Causal variant (true standardised effect = ", BETA_G, ") ---\n", sep = "")
print(data.frame(variant  = colnames(Gmat)[1],
                 neglog10P = round(nl10[1], 3),
                 hat_beta  = round(bhat[1], 4),
                 se_beta   = round(sqrt(vhat[1]), 4)),
      row.names = FALSE)

## calibration on the null variants
nul <- -1L
n_null <- length(nl10) - 1L
p_null <- 10^(-nl10[nul])
lambda <- median(tstat[nul]) / qchisq(0.5, df = 1)

cat("\n--- Calibration on ", n_null, " null variants ---\n", sep = "")
cat("  median -log10(P)          : ", round(median(nl10[nul]), 3),
    "   (expected 0.301)\n", sep = "")
cat("  genomic inflation lambda  : ", round(lambda, 3),
    "   (expected 1.000)\n", sep = "")
cat("  proportion P < 0.05       : ", round(mean(p_null < 0.05), 4),
    "   (expected 0.05)\n", sep = "")
cat("  proportion P < 0.01       : ", round(mean(p_null < 0.01), 4),
    "   (expected 0.01)\n", sep = "")
cat("  mean hat_beta             : ", round(mean(bhat[nul]), 5),
    "   (expected 0)\n", sep = "")

cat("\n--- Timing ---\n")
cat("Whole demo elapsed (s): ",
    round(as.numeric(difftime(Sys.time(), demo_start, units = "secs")), 1),
    "\n", sep = "")
cat("Peak R memory (MB)    : ",
    round(sum(gc()[, "max used"] * c(8, 8) / 1e6)), "\n", sep = "")

cat("\n--- Session ---\n")
cat(R.version.string, "|", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")

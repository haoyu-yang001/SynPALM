## ---------------------------------------------------------------------------
## SynPALM quick-start demo
##
## Self-contained: simulates a small cohort with family structure, so no
## external data are needed. It follows the same call sequence used for the
## UK Biobank analysis in the manuscript:
##
##   build mydf  ->  SynSurrG_ablation_estimate()  ->  score_test_SynSurrG_multiply()
##
## Run with:
##   source(system.file("examples", "quickstart.R", package = "SynPALM"))
##
## The seed is fixed so the printed output matches the "Expected output"
## section of README.md. Do not change it.
## ---------------------------------------------------------------------------

library(SynPALM)
library(Matrix)

set.seed(1)
demo_start <- Sys.time()

## ---------------------------------------------------------------------------
## 1. A sparse GRM with family structure
##
## SynPALM targets cohorts with cryptic relatedness, so the demo uses genuine
## family blocks rather than an identity matrix: 500 families of 4 relatives,
## within-family kinship 0.5 (full siblings).
## ---------------------------------------------------------------------------
n_fam    <- 500
fam_size <- 4
n_all    <- n_fam * fam_size

fam_block <- matrix(0.5, fam_size, fam_size)
diag(fam_block) <- 1
GRM <- Matrix::bdiag(replicate(n_fam, fam_block, simplify = FALSE))
GRM <- methods::as(methods::as(GRM, "CsparseMatrix"), "generalMatrix")

cat("GRM: ", nrow(GRM), " x ", ncol(GRM),
    "  | class: ", class(GRM)[1],
    "  | sparsity: ", round(Matrix::nnzero(GRM) / prod(dim(GRM)), 5), "\n", sep = "")

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
## 3. Genotypes -- 50 variants, the first of which is causal
## ---------------------------------------------------------------------------
n_snp <- 50
maf   <- runif(n_snp, 0.2, 0.4)
Gmat  <- vapply(seq_len(n_snp), function(j) rbinom(n_all, 2, maf[j]), numeric(n_all))
colnames(Gmat) <- paste0("rs_demo_", seq_len(n_snp))

beta_g <- 0.20   # effect size of variant 1, on the standardised scale

## ---------------------------------------------------------------------------
## 4. Target phenotype under the mixed model
##      Y = X beta + g_causal + u + e ,   u ~ N(0, tau * GRM)
## ---------------------------------------------------------------------------
tau   <- 0.4
sigma <- 0.6
L     <- Matrix::t(Matrix::chol(GRM))
u     <- sqrt(tau) * as.vector(L %*% rnorm(n_all))
beta  <- c(0.3, 0.2, rep(0.05, 5))

Y_true <- as.vector(as.matrix(X_all) %*% beta) +
  beta_g * scale(Gmat[, 1])[, 1] +
  u +
  rnorm(n_all, 0, sqrt(sigma))

## ---------------------------------------------------------------------------
## 5. Synthetic / surrogate phenotype
##
## Stands in for the machine-learning prediction (a random forest on clinical
## covariates in the manuscript). SynPALM gains power as its correlation with
## the target rises, and stays valid when that correlation is poor.
## ---------------------------------------------------------------------------
r_ts <- 0.6

int_transform <- function(v) {
  r <- rank(v, na.last = "keep")
  n <- sum(!is.na(v))
  qnorm((r - 0.375) / (n - 2 * 0.375 + 1))
}

S <- int_transform(r_ts * scale(Y_true)[, 1] + sqrt(1 - r_ts^2) * rnorm(n_all))

## ---------------------------------------------------------------------------
## 6. Missingness -- only 20% of the cohort has a measured target phenotype
## ---------------------------------------------------------------------------
obs_frac <- 0.2
obs_idx  <- sort(sample.int(n_all, floor(obs_frac * n_all)))

Y_obs <- rep(NA_real_, n_all)
Y_obs[obs_idx] <- int_transform(Y_true[obs_idx])

cat("Cohort: ", n_all, " individuals | ", length(obs_idx),
    " with an observed target phenotype (", round(100 * obs_frac), "%)\n", sep = "")
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
## Also a sanity check on the block detection: with 500 disjoint families the
## algorithm must recover exactly 500 blocks.
## ---------------------------------------------------------------------------
blocks <- find_blocks_vectorized(GRM)
cat("\nRelatedness blocks found: ", length(blocks),
    " (expected ", n_fam, ")\n", sep = "")
stopifnot(length(blocks) == n_fam)

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
cat("Step 1 elapsed (s): ", round(t_step1[["elapsed"]], 2), "\n", sep = "")

## ---------------------------------------------------------------------------
## 10. Step 2 -- score test across all variants
## ---------------------------------------------------------------------------
cat("\nRunning the SynPALM score test (step 2)...\n")
t_test <- system.time(
  res <- score_test_SynSurrG_multiply(g_matrix = Gmat, step1_pars = SynSurrG_step1)
)
cat("Score test elapsed (s): ", round(t_test[["elapsed"]], 2),
    "  (", n_snp, " variants)\n", sep = "")

## ---------------------------------------------------------------------------
## 11. Results
## ---------------------------------------------------------------------------
cat("\n--- components returned by score_test_SynSurrG_multiply ---\n")
print(names(res))

nl10 <- res[[grep("negative_log10", names(res))[1]]]
bhat <- res[[grep("^hat_beta", names(res))[1]]]

out <- data.frame(
  variant       = colnames(Gmat),
  neglog10P     = round(as.numeric(nl10), 3),
  hat_beta      = round(as.numeric(bhat), 4)
)

cat("\n--- Causal variant ---\n")
print(out[1, ], row.names = FALSE)

cat("\n--- Null variants (2-50) ---\n")
cat("  median -log10(P): ", round(median(out$neglog10P[-1]), 3),
    "   (expected near 0.30)\n", sep = "")
cat("  max    -log10(P): ", round(max(out$neglog10P[-1]), 3), "\n", sep = "")
cat("  proportion with P < 0.05: ",
    round(mean(out$neglog10P[-1] > -log10(0.05)), 3), "\n", sep = "")

cat("\n--- Timing ---\n")
cat("Whole demo elapsed (s): ",
    round(as.numeric(difftime(Sys.time(), demo_start, units = "secs")), 2),
    "\n", sep = "")

cat("\n--- Session ---\n")
cat(R.version.string, "|", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")

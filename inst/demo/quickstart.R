## ---------------------------------------------------------------------------
## SynPALM quick-start demo
##
## Self-contained: simulates a small dataset, so no external data are needed.
## Run with:
##   source(system.file("demo", "quickstart.R", package = "SynPALM"))
##
## The seed below is fixed so that the printed output matches the "Expected
## output" section of README.md exactly. Do not change it.
## ---------------------------------------------------------------------------

library(SynPALM)
library(Matrix)

set.seed(1)

demo_start <- Sys.time()

## 1. Simulation parameters -------------------------------------------------
n_obs        <- 1000   # samples with an observed primary phenotype
missing_rate <- 0.5    # 50% missingness in the target phenotype
maf          <- 0.3    # minor allele frequency
rho          <- 0.5    # within-sibling correlation

## 2. Genotypes -------------------------------------------------------------
G_matrix <- generate_genotypes_rho(
  n_obs  = n_obs,
  miss   = missing_rate,
  n_snps = 1,
  maf    = maf,
  rho    = rho
)

## 3. Genetic relatedness matrix -------------------------------------------
n_all <- length(G_matrix)
GRM   <- as(diag(n_all), "dgCMatrix")

## 4. Step-1 quantities (mocked for this demo) ------------------------------
pre_vars <- list(
  chol_Sigma22              = diag(n_all),
  L_cond                    = diag(n_obs),
  Sigma12_Sigma22inv        = matrix(0.5, n_obs, n_all),
  Sigma12_oracle_Sigma22inv = matrix(0.5, n_all, n_all),
  L_cond_oracle             = diag(n_all)
)

## 5. Simulate phenotypes --------------------------------------------------
sim_data <- DGP_GRMr_gmodel_step3_rho(
  n_obs           = n_obs,
  miss            = missing_rate,
  pve_g           = 0.005,
  Genotype_each   = G_matrix,
  pre_vars_step2  = pre_vars,
  GRM             = GRM
)

## 6. Syn-PALM score test --------------------------------------------------
timing  <- system.time(
  results <- score_test_SynSurrG_power_sim_one_g_to_one_Y(sim_data)
)

## 7. Report ---------------------------------------------------------------
cat("\n--- Syn-PALM demo output ---\n")
print(results)

cat("\n--- Timing ---\n")
cat("Score test elapsed (s): ", round(timing[["elapsed"]], 2), "\n", sep = "")
cat("Whole demo elapsed (s): ",
    round(as.numeric(difftime(Sys.time(), demo_start, units = "secs")), 2),
    "\n", sep = "")

cat("\n--- Session ---\n")
cat(R.version.string, "|", Sys.info()[["sysname"]], Sys.info()[["release"]], "\n")

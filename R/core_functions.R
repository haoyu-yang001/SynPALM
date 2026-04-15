## functions ########
library(lattice)
library(data.table)
library(dplyr)
library(Matrix)
library(GMMAT)
library(ranger)
library(parallel)
library(withr)



#' Syn-PALM Score Test for Multiple SNPs
#'
#' This function implements the core score test of the Syn-PALM (Synthetic Phenotype
#' Assisted Linear Mixed models) framework. It jointly analyzes observed proteomic
#' measurements and synthetic data for a batch of SNPs, accounting for sample
#' relatedness using pre-computed parameters from Step 1.
#'
#' @param g_matrix A numeric genotype matrix (n x m) where n is the total number
#'   of samples and m is the number of SNPs.
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model fitting.
#'   Expected elements include \code{inv_Sigma22}, \code{obs_protein_index}, \code{V11},
#'   \code{hatY}, \code{Y}, \code{X_obs}, \code{Att}, and others.
#'
#' @return A list containing the following components for each SNP:
#' \itemize{
#'   \item{T_score_SynSurrG}: The Syn-PALM score test statistics.
#'   \item{negative_log10_pval_SynSurrG}: P-values on the -log10 scale.
#'   \item{hat_beta_SynSurrG}: Estimated effect sizes (beta).
#'   \item{var_hat_beta_SynSurrG}: Estimated variances of the effect sizes.
#' }
#' @export
#' @import Matrix
score_test_SynSurrG_multiply <- function(g_matrix, step1_pars) {

  # Attach parameters from Step 1 to the current environment
  list2env(step1_pars, envir = environment())

  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)

  snpindex <- 1:ncol(g_matrix)

  G_all <- g_matrix
  G_obs <- as.matrix(G_all[obs_protein_index,])

  # Efficient matrix operations for Score and Variance calculation
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X

  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)

  bt1 <- Sigma12_Sigma22inv %*% G_all
  Ah <- t(temA) %*% t(Sigma12)
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)

  # Internal function to calculate statistics for a single SNP
  compute_snp <- function(snp_d) {

    ## 1. Estimate alpha (Nuisance parameters for synthetic data)
    # Construct block matrix A
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)

    # Construct the right-hand side for alpha estimation
    Ar0 <- rbind(t(temA[,snp_d]), t(Att)) # Z %*% Sigma_22^{-1}

    # Solve for alpha estimates
    alpha <- solve(A_mat, Ar0 %*% hatY)

    ## 2. Estimate beta
    b1 <- cbind(bt1[,snp_d], bt2)
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)

    ## 3. Calculate the score (S_syn)
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))

    ## 4. Calculate the variance of the score (VU)
    # VU involves complex combinations of variance-covariance components
    Arr <- rbind(Ah[snp_d,], Atb)
    Ar <- solve(A_mat, Arr)

    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    V2 <- VV2 %*% VV1[snp_d,]

    VU <- as.numeric(V1[snp_d] + V2)

    # Effect size and its variance estimation
    inv_beta <- VV1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)

    list(S_syn = S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }

  # Apply computation over all SNPs in the batch
  results <- lapply(seq_along(snpindex), compute_snp)

  # Extract and consolidate results
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")

  # Compute final scores and p-values
  surr_result <- compute_score(S_syn, VU_syn)

  return(list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  ))
}

#' Score Test for Multiple SNPs using the SynSurr Method
#'
#' This function performs the SynSurr (Synthetic Surrogate) score test for a batch of SNPs.
#' Unlike Syn-PALM which jointly analyzes observed and synthetic data, SynSurr
#' focuses on the analysis utilizing synthetic phenotypes as a surrogate,
#' while accounting for sample relatedness and population structure.
#'
#' @param g_matrix A numeric genotype matrix (n_total x m).
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model fitting
#'   (e.g., \code{inv_Sigma22}, \code{V11}, \code{hatY}, etc.).
#' @param independent_indices An integer vector of row indices used to subset
#'   \code{g_matrix} to ensure alignment with the analysis sample.
#'
#' @return A list containing the following components for each SNP:
#' \itemize{
#'   \item{T_score_SynSurr}: The SynSurr score test statistics.
#'   \item{negative_log10_pval_SynSurr}: P-values on the -log10 scale.
#'   \item{hat_beta_SynSurr}: Estimated effect sizes (beta).
#'   \item{var_hat_beta_SynSurr}: Estimated variances of the effect sizes.
#' }
#' @export
#' @import Matrix
score_test_SynSurr_multiply <- function(g_matrix, step1_pars, independent_indices) {

  # Load Step 1 pre-computed parameters into the current function environment
  list2env(step1_pars, envir = environment())

  # Subset genotype matrix to match the independent sample indices
  g_matrix <- as.matrix(g_matrix[independent_indices,])

  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)

  snpindex <- 1:ncol(g_matrix)

  G_all <- g_matrix
  G_obs <- as.matrix(G_all[obs_protein_index,])

  # Efficient matrix computations for Score and Variance components
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X

  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)

  bt1 <- Sigma12_Sigma22inv %*% G_all

  Ah <- t(temA) %*% t(Sigma12)

  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)

  V1 <- rowSums((VV1 %*% Sigma11) * VV1)

  # Define internal function for single SNP computation
  compute_snp <- function(snp_d) {

    ## 1. Estimate alpha (Nuisance parameters for the synthetic component)
    # Construct block matrix A
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)

    # Construct right-hand side for alpha estimation
    Ar0 <- rbind(t(temA[,snp_d]), t(Att)) # Z %*% Sigma_22^{-1}

    # Solve for alpha estimates
    alpha <- solve(A_mat, Ar0 %*% hatY)

    ## 2. Estimate beta
    b1 <- cbind(bt1[,snp_d], bt2)
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)

    ## 3. Calculate the score (S_syn)
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))

    ## 4. Calculate the variance of the score (VU)
    Arr <- rbind(Ah[snp_d,], Atb)
    Ar <- solve(A_mat, Arr)

    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar

    V2 <- VV2 %*% VV1[snp_d,]

    VU <- as.numeric(V1[snp_d] + V2)

    # Calculate effect size estimate and its variance
    inv_beta <- VV1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)

    list(S_syn = S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }

  # Apply computation across all SNPs
  results <- lapply(seq_along(snpindex), compute_snp)

  # Consolidate results into vectors/matrices
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")

  # Compute final score statistics and p-values
  surr_result <- compute_score(S_syn, VU_syn)

  return(list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  ))
}

#' Score Test for Multiple SNPs using the Oracle Method
#'
#' This function performs the "Oracle" score test for a batch of SNPs. The Oracle
#' method serves as a benchmark by assuming all phenotypes are observed for all
#' samples. It calculates the score statistics and effect size estimates within
#' a Linear Mixed Model (LMM) framework, accounting for sample relatedness.
#'
#' @param g_matrix A numeric genotype matrix (n x m) for all samples.
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model fitting,
#'   specifically those derived under the assumption of fully observed data (e.g.,
#'   \code{invS11_res}, \code{inv_Sigma11}, \code{SX}, \code{XtSX_inv}).
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_OracleG}: The Oracle score test statistics.
#'   \item{negative_log10_pval_OracleG}: P-values on the -log10 scale.
#'   \item{hat_beta_OracleG}: Estimated effect sizes (beta) for the Oracle model.
#'   \item{var_hat_beta_OracleG}: Estimated variances of the effect sizes.
#' }
#' @export
score_test_OracleG_multiply <- function(g_matrix, step1_pars) {

  # Load pre-computed parameters into the current environment
  list2env(step1_pars, envir = environment())

  ## Calculate Oracle test statistics
  # SU: Score vector (numerator)
  SU <- as.numeric(crossprod(g_matrix, invS11_res))

  # A: Intermediate projection of genotypes onto the covariate space
  A <- crossprod(SX, g_matrix)

  # VU: Variance of the score vector (denominator)
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))

  # Compute p-values based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Oracle method
  names(results) <- c("T_score_OracleG", "negative_log10_pval_OracleG")

  # Calculate effect size (beta) and its variance
  # P1G represents the projection adjustment
  P1G <- X %*% (XtSX_inv %*% A)

  # inv_beta: The information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix))) -
    as.numeric(colSums(g_matrix * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_OracleG <- hat_beta
  results$var_hat_beta_OracleG <- var_hat_beta

  return(results)
}


#' Score Test for Multiple SNPs using the Oracle Method (Independent Subset)
#'
#' This function performs the "Oracle" score test for a batch of SNPs using only
#' a subset of independent samples. It serves as a benchmark to evaluate the
#' performance gain of Syn-PALM over standard analysis on independent individuals,
#' assuming their phenotypes are fully observed.
#'
#' @param g_matrix A numeric genotype matrix (n_total x m).
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model
#'   fitting for independent samples (e.g., \code{invS11_res}, \code{inv_Sigma11},
#'   \code{SX}, \code{XtSX_inv}).
#' @param independent_indices An integer vector of row indices identifying the
#'   independent subset of samples to be included in the analysis.
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_Oracle}: The Oracle score test statistics for independent samples.
#'   \item{negative_log10_pval_Oracle}: P-values on the -log10 scale.
#'   \item{hat_beta_Oracle}: Estimated effect sizes (beta) for the Oracle model.
#'   \item{var_hat_beta_Oracle}: Estimated variances of the effect sizes.
#' }
#' @export
score_test_Oracle_multiply <- function(g_matrix, step1_pars, independent_indices) {

  # Load pre-computed Step 1 parameters into the current environment
  list2env(step1_pars, envir = environment())

  # Subset genotype matrix to include only independent individuals
  g_matrix <- as.matrix(g_matrix[independent_indices,])

  ## Calculate Oracle test statistics for independent samples
  # SU: Score vector (numerator)
  SU <- as.numeric(crossprod(g_matrix, invS11_res))

  # A: Projection of genotypes onto the covariate space
  A <- crossprod(SX, g_matrix)

  # VU: Variance of the score vector (denominator)
  # Calculated as the difference between total variance and variance explained by covariates
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))

  # Compute p-values based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Oracle (Independent) method
  names(results) <- c("T_score_Oracle", "negative_log10_pval_Oracle")

  # Calculate effect size (beta) and its variance
  # P1G represents the projection adjustment for covariates
  P1G <- X %*% (XtSX_inv %*% A)

  # inv_beta: The denominator for the beta estimator (information content)
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix))) -
    as.numeric(colSums(g_matrix * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_Oracle <- hat_beta
  results$var_hat_beta_Oracle <- var_hat_beta

  return(results)
}

#' Score Test for Multiple SNPs using the Observed-only Method (with Relatedness)
#'
#' This function performs a score test for a batch of SNPs using only the samples
#' with observed proteomic measurements. It accounts for cryptic relatedness
#' and population structure among the observed individuals using a Linear
#' Mixed Model (LMM) framework.
#'
#' @param g_matrix A numeric genotype matrix (n_total x m).
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model
#'   fitting, specifically restricted to or derived from the observed individuals
#'   (e.g., \code{obs_protein_index}, \code{invS11_res}, \code{inv_Sigma11},
#'   \code{SX}, \code{XtSX_inv}).
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_ObsG}: The score test statistics based on observed data.
#'   \item{negative_log10_pval_ObsG}: P-values on the -log10 scale.
#'   \item{hat_beta_ObsG}: Estimated effect sizes (beta) using observed data.
#'   \item{var_hat_beta_ObsG}: Estimated variances of the effect sizes.
#' }
#' @export
score_test_ObsG_multiply <- function(g_matrix, step1_pars) {

  # Load Step 1 parameters into the local environment
  list2env(step1_pars, envir = environment())

  # Subset genotype matrix to include only samples with observed protein data
  g_matrix_sub <- g_matrix[obs_protein_index, ]

  ## Calculate test statistics for observed data
  # SU: Score vector (numerator)
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))

  # A: Projection of subsetted genotypes onto the covariate space
  A <- crossprod(SX, g_matrix_sub)

  # VU: Variance of the score vector (denominator) accounting for relatedness
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) -
                     colSums(A * (XtSX_inv %*% A)))

  # Compute p-values based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Observed-only (with relatedness) method
  names(results) <- c("T_score_ObsG", "negative_log10_pval_ObsG")

  # Calculate effect size (beta) and its variance
  # P1G: Projection adjustment for the subsetted data
  P1G <- X %*% (XtSX_inv %*% A)

  # inv_beta: Information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) -
    as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta

  return(results)
}

#' Score Test for Multiple SNPs using the Observed-only Method (Independent Subset)
#'
#' This function performs a score test for a batch of SNPs using only the
#' independent samples that have observed proteomic measurements. It does not
#' account for cryptic relatedness (standard regression framework), serving as
#' a baseline to evaluate the impact of including synthetic data and
#' accounting for sample relatedness.
#'
#' @param g_matrix A numeric genotype matrix (n_total x m).
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model
#'   fitting for independent observed samples (e.g., \code{invS11_res},
#'   \code{inv_Sigma11}, \code{SX}, \code{XtSX_inv}).
#' @param independent_indices An integer vector of row indices identifying the
#'   subset of independent samples in the study.
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_Obs}: The score test statistics based on independent observed data.
#'   \item{negative_log10_pval_Obs}: P-values on the -log10 scale.
#'   \item{hat_beta_Obs}: Estimated effect sizes (beta) for the observed independent samples.
#'   \item{var_hat_beta_Obs}: Estimated variances of the effect sizes.
#' }
#' @export
score_test_Obs_multiply <- function(g_matrix, step1_pars, independent_indices) {

  # Load Step 1 parameters into the local environment
  list2env(step1_pars, envir = environment())

  # First, subset genotype matrix to include only independent individuals
  g_matrix_all_sub <- as.matrix(g_matrix[independent_indices,])

  # Second, subset further to include only those with observed protein data
  g_matrix_sub <- g_matrix_all_sub[obs_protein_index,]

  ## Calculate test statistics for observed independent data
  # SU: Score vector (numerator)
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))

  # A: Projection of subsetted genotypes onto the covariate space
  A <- crossprod(SX, g_matrix_sub)

  # VU: Variance of the score vector (denominator)
  # In this context, inv_Sigma11 usually represents a diagonal matrix (no relatedness)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) -
                     colSums(A * (XtSX_inv %*% A)))

  # Compute p-values based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Observed-only (Independent) method
  names(results) <- c("T_score_Obs", "negative_log10_pval_Obs")

  # Calculate effect size (beta) and its variance
  # P1G: Projection adjustment for the subsetted data
  P1G <- X %*% (XtSX_inv %*% A)

  # inv_beta: Information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) -
    as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results
  results$hat_beta_Obs <- hat_beta
  results$var_hat_beta_Obs <- var_hat_beta

  return(results)
}

#' Syn-PALM Score Test for a Single SNP
#'
#' This function implements the Syn-PALM (Synthetic Phenotype Assisted Linear
#' Mixed models) score test specifically optimized for a single genetic variant.
#' It computes the score statistic, p-value, and effect size by jointly
#' analyzing observed proteomic data and complete synthetic data.
#'
#' @param G_all A numeric vector representing the genotypes for a single SNP
#'   across all $n$ samples.
#' @param step1_pars A list containing pre-computed model parameters and
#'   variance-covariance components from the Step 1 fitting process.
#'
#' @return A list containing the following components for the SNP:
#' \itemize{
#'   \item{T_score_SynSurrG}: The Syn-PALM score test statistic.
#'   \item{negative_log10_pval_SynSurrG}: P-value on the -log10 scale.
#'   \item{hat_beta_SynSurrG}: Estimated effect size (beta).
#'   \item{var_hat_beta_SynSurrG}: Estimated variance of the effect size.
#' }
#' @export
#' @import Matrix
score_test_SynSurrG_single <- function(G_all, step1_pars) {

  # Load Step 1 parameters into the local function environment
  list2env(step1_pars, envir = environment())

  n_obs <- length(obs_protein_index)
  n <- length(G_all)

  # Extract genotypes for samples with observed protein measurements
  G_obs <- G_all[obs_protein_index]

  # Pre-calculate components for the Score and Variance computation
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X

  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)

  bt1 <- Sigma12_Sigma22inv %*% G_all

  Ah <- t(temA) %*% t(Sigma12)

  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)

  V1 <- rowSums((VV1 %*% Sigma11) * VV1)

  # Internal function to compute statistics for the variant
  compute_snp <- function(snp_d) {

    ## 1. Estimate alpha (Nuisance parameters for the synthetic component)
    # Construct block matrix A for the joint model
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)

    # Construct the right-hand side for alpha estimation
    Ar0 <- rbind(t(temA[,snp_d]), t(Att)) # Z %*% Sigma_22^{-1}

    # Solve for alpha estimates
    alpha <- solve(A_mat, Ar0 %*% hatY)

    ## 2. Estimate beta (Effect size)
    b1 <- cbind(bt1[,snp_d], bt2)
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)

    ## 3. Calculate the score (S_syn)
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))

    ## 4. Calculate the variance of the score (VU)
    Arr <- rbind(Ah[snp_d,], Atb)
    Ar <- solve(A_mat, Arr)

    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar

    V2 <- VV2 %*% VV1[snp_d,]

    VU <- as.numeric(V1[snp_d] + V2)

    # Calculate final effect size and variance estimates
    inv_beta <- VV1 %*% G_obs
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)

    list(S_syn = S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }

  # Execute computation
  results <- lapply(seq_along(1), compute_snp)

  # Consolidate results
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")

  # Compute final test statistics and p-values
  surr_result <- compute_score(S_syn, VU_syn)

  return(list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  ))
}

#' SynSurr Score Test for a Single SNP (Independent Subset)
#'
#' This function performs the SynSurr (Synthetic Surrogate) score test for a
#' single genetic variant. It utilizes independent sample indices to subset the
#' data and calculates the score statistic, p-value, and effect size based on
#' synthetic surrogate phenotypes.
#'
#' @param G_all A numeric vector containing the genotypes for a single SNP
#'   across all $n$ samples.
#' @param step1_pars A list of pre-calculated model parameters and
#'   variance-covariance components from the Step 1 fitting process.
#' @param independent_indices An integer vector of row indices used to subset
#'   \code{G_all} to include only independent individuals in the analysis.
#'
#' @return A list containing the following components for the SNP:
#' \itemize{
#'   \item{T_score_SynSurr}: The SynSurr score test statistic.
#'   \item{negative_log10_pval_SynSurr}: P-value on the -log10 scale.
#'   \item{hat_beta_SynSurr}: Estimated effect size (beta).
#'   \item{var_hat_beta_SynSurr}: Estimated variance of the effect size.
#' }
#' @export
#' @import Matrix
score_test_SynSurr_single <- function(G_all, step1_pars, independent_indices) {

  # Load Step 1 parameters into the local environment
  list2env(step1_pars, envir = environment())

  # Subset genotype vector for independent individuals
  G_all <- G_all[independent_indices]

  n_obs <- length(obs_protein_index)
  n <- length(G_all)

  # Extract genotypes for individuals with observed protein data
  G_obs <- G_all[obs_protein_index]

  # Pre-compute matrix components for score and variance calculation efficiency
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X

  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)

  bt1 <- Sigma12_Sigma22inv %*% G_all

  Ah <- t(temA) %*% t(Sigma12)

  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)

  V1 <- rowSums((VV1 %*% Sigma11) * VV1)

  # Internal function for calculating statistics for the single variant
  compute_snp <- function(snp_d) {

    ## 1. Estimate alpha (Nuisance parameters for synthetic phenotypes)
    # Construct block matrix A
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)

    # Construct right-hand side for alpha estimation
    Ar0 <- rbind(t(temA[,snp_d]), t(Att)) # Z %*% Sigma_22^{-1}

    # Solve for alpha estimates
    alpha <- solve(A_mat, Ar0 %*% hatY)

    ## 2. Estimate beta (Effect size)
    b1 <- cbind(bt1[,snp_d], bt2)
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)

    ## 3. Calculate the score (S_syn)
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))

    ## 4. Calculate the variance of the score (VU)
    Arr <- rbind(Ah[snp_d,], Atb)
    Ar <- solve(A_mat, Arr)

    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar

    V2 <- VV2 %*% VV1[snp_d,]

    VU <- as.numeric(V1[snp_d] + V2)

    # Calculate final effect size and variance estimates
    inv_beta <- VV1 %*% G_obs
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)

    list(S_syn = S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }

  # Execute the computation (using lapply for consistency with batch functions)
  results <- lapply(seq_along(1), compute_snp)

  # Consolidate results into vectors/scalars
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")

  # Compute final score test results and p-values
  surr_result <- compute_score(S_syn, VU_syn)

  return(list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  ))
}

#' Oracle Score Test for a Single SNP (with Relatedness)
#'
#' This function performs the Oracle score test for a single genetic variant
#' across all samples, assuming phenotypes are fully observed. It accounts for
#' cryptic relatedness and population structure using a Linear Mixed Model (LMM)
#' framework, providing an ideal benchmark for statistical power.
#'
#' @param G_all A numeric vector containing the genotypes for a single SNP
#'   across all $n$ samples.
#' @param step1_pars A list of pre-calculated parameters from the Step 1
#'   LMM fitting (e.g., \code{invS11_res}, \code{inv_Sigma11}, \code{SX},
#'   \code{XtSX_inv}, \code{X}).
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_OracleG}: The Oracle score test statistic for the SNP.
#'   \item{negative_log10_pval_OracleG}: P-value on the -log10 scale.
#'   \item{hat_beta_OracleG}: Estimated effect size (beta) under the Oracle model.
#'   \item{var_hat_beta_OracleG}: Estimated variance of the effect size.
#' }
#' @export
score_test_OracleG_single <- function(G_all, step1_pars) {

  # Load pre-computed parameters into the local environment
  list2env(step1_pars, envir = environment())

  ## Calculate Oracle test statistics for the single variant
  # SU: Score statistic numerator
  SU <- as.numeric(crossprod(G_all, invS11_res))

  # A: Intermediate projection of genotypes onto covariate space
  A <- crossprod(SX, G_all)

  # VU: Score statistic denominator (variance) accounting for relatedness
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) -
                     colSums(A * (XtSX_inv %*% A)))

  # Compute p-value using the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs to reflect the Oracle method
  names(results) <- c("T_score_OracleG", "negative_log10_pval_OracleG")

  # Calculate effect size (beta) and its variance
  # P1G: Projection adjustment representing the part of G explained by X
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all)

  # inv_beta: Information content for the single variant
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all))) -
    as.numeric(colSums(G_all * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_OracleG <- hat_beta
  results$var_hat_beta_OracleG <- var_hat_beta

  return(results)
}

#' Oracle Score Test for a Single SNP (Independent Subset)
#'
#' This function performs the Oracle score test for a single genetic variant using
#' a subset of independent samples. It assumes phenotypes are fully observed for
#' these individuals and does not account for cryptic relatedness, providing a
#' baseline benchmark for standard regression-based analysis.
#'
#' @param G_all A numeric vector containing the genotypes for a single SNP
#'   across all $n$ samples.
#' @param step1_pars A list of pre-calculated parameters from the Step 1 model
#'   fitting for independent samples (e.g., \code{invS11_res}, \code{inv_Sigma11},
#'   \code{SX}, \code{XtSX_inv}, \code{X}).
#' @param independent_indices An integer vector of row indices identifying the
#'   independent subset of samples to be included in the analysis.
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_Oracle}: The Oracle score test statistic for the independent subset.
#'   \item{negative_log10_pval_Oracle}: P-value on the -log10 scale.
#'   \item{hat_beta_Oracle}: Estimated effect size (beta) for the Oracle model.
#'   \item{var_hat_beta_Oracle}: Estimated variance of the effect size.
#' }
#' @export
score_test_Oracle_single <- function(G_all, step1_pars, independent_indices) {

  # Load pre-computed parameters into the local environment
  list2env(step1_pars, envir = environment())

  # Subset genotype vector to include only independent individuals
  G_all <- G_all[independent_indices]

  ## Calculate Oracle test statistics for the independent variant
  # SU: Score statistic numerator
  SU <- as.numeric(crossprod(G_all, invS11_res))

  # A: Intermediate projection of genotypes onto covariate space
  A <- crossprod(SX, G_all)

  # VU: Score statistic denominator (variance)
  # Calculated as the difference between total variance and variance explained by covariates
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))

  # Compute p-value based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Oracle (Independent) method
  names(results) <- c("T_score_Oracle", "negative_log10_pval_Oracle")

  # Calculate effect size (beta) and its variance
  # P1G represents the projection of G onto the covariate space X
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all)

  # inv_beta: Information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all))) -
    as.numeric(colSums(G_all * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_Oracle <- hat_beta
  results$var_hat_beta_Oracle <- var_hat_beta

  return(results)
}

#' Observed-only Score Test for a Single SNP (with Relatedness)
#'
#' This function performs a score test for a single genetic variant using only the
#' samples with observed proteomic measurements. It accounts for cryptic
#' relatedness and population structure among the observed individuals using a
#' Linear Mixed Model (LMM) framework.
#'
#' @param G_all A numeric vector containing the genotypes for a single SNP
#'   across all samples.
#' @param step1_pars A list of pre-calculated parameters from the Step 1
#'   LMM fitting, specifically restricted to or derived from the observed
#'   individuals (e.g., \code{obs_protein_index}, \code{invS11_res},
#'   \code{inv_Sigma11}, \code{SX}, \code{XtSX_inv}, \code{X}).
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_ObsG}: The score test statistic based on observed data for the SNP.
#'   \item{negative_log10_pval_ObsG}: P-value on the -log10 scale.
#'   \item{hat_beta_ObsG}: Estimated effect size (beta) using observed data.
#'   \item{var_hat_beta_ObsG}: Estimated variance of the effect size.
#' }
#' @export
score_test_ObsG_single <- function(G_all, step1_pars) {

  # Load pre-computed Step 1 parameters into the current environment
  list2env(step1_pars, envir = environment())

  # Subset genotype vector to include only samples with observed protein data
  G_all_sub <- G_all[obs_protein_index]

  ## Calculate test statistics for observed data
  # SU: Score statistic numerator
  SU <- as.numeric(crossprod(G_all_sub, invS11_res))

  # A: Intermediate projection of subsetted genotypes onto covariate space
  A <- crossprod(SX, G_all_sub)

  # VU: Score statistic denominator (variance) accounting for relatedness
  VU <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub)) -
                     colSums(A * (XtSX_inv %*% A)))

  # Compute p-value based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Observed-only (with relatedness) method
  names(results) <- c("T_score_ObsG", "negative_log10_pval_ObsG")

  # Calculate effect size (beta) and its variance
  # P1G: Projection adjustment representing the part of G_sub explained by X
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all_sub)

  # inv_beta: Information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub))) -
    as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta

  return(results)
}

#' Observed-only Score Test for a Single SNP (Independent Subset)
#'
#' This function performs a score test for a single genetic variant using only
#' the independent samples that have observed proteomic measurements. This
#' method does not account for cryptic relatedness (standard regression
#' framework) and serves as a baseline for comparing the efficiency of the
#' Syn-PALM approach.
#'
#' @param G_all A numeric vector containing the genotypes for a single SNP
#'   across all samples.
#' @param step1_pars A list of pre-calculated parameters from the Step 1
#'   fitting process for independent observed samples (e.g., \code{invS11_res},
#'   \code{inv_Sigma11}, \code{SX}, \code{XtSX_inv}, \code{X}).
#' @param independent_indices An integer vector of row indices identifying the
#'   subset of independent samples in the study.
#'
#' @return A list containing:
#' \itemize{
#'   \item{T_score_Obs}: The score test statistic based on independent observed data.
#'   \item{negative_log10_pval_Obs}: P-value on the -log10 scale.
#'   \item{hat_beta_Obs}: Estimated effect size (beta) for the observed independent samples.
#'   \item{var_hat_beta_Obs}: Estimated variance of the effect size.
#' }
#' @export
score_test_Obs_single <- function(G_all, step1_pars, independent_indices) {

  # Load Step 1 parameters into the local environment
  list2env(step1_pars, envir = environment())

  # First, subset genotype vector to include only independent individuals
  G_all_sub <- G_all[independent_indices]

  # Second, subset further to include only those with observed protein data
  G_sub <- G_all_sub[obs_protein_index]

  ## Calculate test statistics for observed independent data
  # SU: Score statistic numerator
  SU <- as.numeric(crossprod(G_sub, invS11_res))

  # A: Intermediate projection of subsetted genotypes onto covariate space
  A <- crossprod(SX, G_sub)

  # VU: Score statistic denominator (variance)
  # In this context, inv_Sigma11 usually represents a diagonal matrix
  VU <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub)) -
                     colSums(A * (XtSX_inv %*% A)))

  # Compute p-value based on the score and variance
  results <- compute_score(SU, VU)

  # Rename outputs for the Observed-only (Independent) method
  names(results) <- c("T_score_Obs", "negative_log10_pval_Obs")

  # Calculate effect size (beta) and its variance
  # P1G: Projection adjustment representing the part of G_sub explained by X
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_sub)

  # inv_beta: Information content/denominator for the beta estimator
  inv_beta <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub))) -
    as.numeric(colSums(G_sub * (inv_Sigma11 %*% P1G)))

  hat_beta <- SU / (inv_beta)
  var_hat_beta <- VU / (inv_beta)^2

  # Append beta estimates to the results list
  results$hat_beta_Obs <- hat_beta
  results$var_hat_beta_Obs <- var_hat_beta

  return(results)
}


#' Step 1 Parameter Estimation for Syn-PALM (with Relatedness)
#'
#' This function performs the first step of the Syn-PALM pipeline for GWAS simulation.
#' It estimates variance components for both target and surrogate phenotypes,
#' constructs the joint covariance matrices (Sigma11, Sigma12, Sigma22), and
#' pre-calculates essential components for the SynSurr score test using
#' Generalized Least Squares (GLS) logic.
#'
#' @param mydf A list containing the simulated data (usually the output of \code{DGP}):
#'   \itemize{
#'     \item{\code{X_all}}: Covariates for all samples.
#'     \item{\code{Y}}: True (Oracle) phenotypes (not used for estimation).
#'     \item{\code{Y_obs}}: Partially observed target phenotypes.
#'     \item{\code{S}}: Synthetic/Surrogate phenotypes.
#'     \item{\code{GRM}}: Genetic Relatedness Matrix.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) required for
#'   the single-SNP score test, including:
#' \itemize{
#'   \item{\code{V11}}: The weight matrix for the score numerator.
#'   \item{\code{B_mat}}: Projection matrix for beta estimation.
#'   \item{\code{Sigma12_Sigma22inv}}: The projection of surrogate data onto the target space.
#'   \item{... and other pre-computed matrix components.}
#' }
#' @export
#' @import Matrix
SynSurrG_typeIerror_step1 <- function(mydf) {

  # Organize input data
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  GRM <- mydf$GRM
  n <- nrow(test_df)

  # Define indices for observed protein data
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index, 1:n]
  n_obs <- length(obs_protein_index)

  # Prepare data for variance component estimation
  reml_data_SynSurrG <- na.omit(data.frame(y = test_df$Y_obs, haty = test_df$yhat, x = test_df$X))
  lm_data_SynSurrG <- data.frame(y = test_df$Y_obs, haty = test_df$yhat, x = test_df$X)

  # Fit initial null models to obtain residuals
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurrG)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurrG)

  # Define off-diagonal GRMs for Method-of-Moments estimation
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0

  ## Method-of-Moments Variance Component Estimation
  # Target Protein (T)
  hat_tau_T2 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals) / sum((GRM_o)^2)
  hat_sigma_T2 <- (sum(model_lm1$residuals^2) - hat_tau_T2 * sum(diag(GRM_obs_obs))) / n_obs

  # Surrogate Protein (S)
  hat_tau_S2 <- as.numeric(t(model_lm2$residuals) %*% GRM_oall %*% model_lm2$residuals) / sum((GRM_oall)^2)
  hat_sigma_S2 <- (sum(model_lm2$residuals^2) - hat_tau_S2 * sum(diag(GRM))) / n

  # Cross-trait (TS)
  hat_tau_TS <- as.numeric(t(model_lm2$residuals[obs_protein_index]) %*% GRM_o %*% model_lm1$residuals) / sum((GRM_o)^2)
  hat_sigma_TS <- (sum(model_lm2$residuals[obs_protein_index] * model_lm1$residuals) -
                     hat_tau_TS * sum(diag(GRM_obs_obs))) / n_obs

  ## Construct Covariance Matrices
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)

  # Construct Sigma12 (Cross-covariance)
  I_matrix <- sparseMatrix(i = 1:min(n_obs, n), j = 1:min(n_obs, n),
                           x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n))
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix

  # Invert matrices using block-wise optimization
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)
  inv_Sigma22 <- matrix_inv_block(wait_matrix = Sigma22)

  ## Pre-calculate SynSurr components
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Calculate A based on Schur Complement logic
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- matrix_inv_Amatrix(Amatrix = A)

  # Define observed and synthetic phenotypes for Step 2
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat

  # Define covariate matrices (including intercept)
  X_all <- as.matrix(cbind(intercept = 1, X = test_df$X))
  X_obs <- X_all[obs_protein_index, ]

  # Pre-compute SNP-independent algebraic terms
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))

  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))

  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all

  # Export all pre-calculated parameters
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att",
              "A22", "Atb", "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
              "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12")

  return(mget(params, envir = environment()))
}

#' Step 1 Parameter Estimation for Syn-PALM Ablation Study
#'
#' This function performs the Step 1 estimation process specifically for ablation
#' studies within the Syn-PALM framework. It supports multiple covariates by
#' dynamically detecting columns starting with "X." and ensures numerical stability
#' in variance component estimation through a lower-bound constraint on residual variance.
#'
#' @param mydf A list containing the simulated or real data:
#'   \itemize{
#'     \item{\code{X_all}}: A data frame or matrix of covariates (columns prefixed with "X.").
#'     \item{\code{Y}}: Oracle phenotypes.
#'     \item{\code{Y_obs}}: Partially observed target phenotypes.
#'     \item{\code{S}}: Synthetic/Surrogate phenotypes.
#'     \item{\code{GRM}}: Genetic Relatedness Matrix.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) optimized for
#'   multi-covariate models, including weight matrices (\code{V11}), projection
#'   matrices (\code{B_mat}), and joint covariance components.
#' @export
#' @import Matrix
SynSurrG_ablation_estimate <- function(mydf) {

  # Organize input data and identify covariate columns
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  n <- nrow(test_df)
  GRM <- mydf$GRM

  # Define indices for observed individuals
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index, 1:n]
  n_obs <- length(obs_protein_index)

  # Prepare datasets for model fitting (handling multi-covariates via regex)
  reml_data_SynSurrG <- na.omit(data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    test_df[, grepl("^X\\.", names(test_df))]
  ))

  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    test_df[, grepl("^X\\.", names(test_df))]
  )

  # Dynamically construct model formulas for Y and Surrogate S
  preds <- grep("^x\\.X\\.", names(reml_data_SynSurrG), value = TRUE)
  fml_y <- reformulate(preds, response = "y")
  model_lm1 <- lm(fml_y, data = reml_data_SynSurrG)

  fml_haty <- reformulate(preds, response = "haty")
  model_lm2 <- lm(fml_haty, data = lm_data_SynSurrG)

  # Variance Component Estimation using Method-of-Moments
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0

  # Target variance (tau_T and sigma_T) - Includes a 0.01 floor for stability
  hat_tau_T2 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals) / sum((GRM_o)^2)
  hat_sigma_T2 <- max((sum(model_lm1$residuals^2) - hat_tau_T2 * sum(diag(GRM_obs_obs))) / n_obs, 0.01)

  # Surrogate variance (tau_S and sigma_S)
  hat_tau_S2 <- as.numeric(t(model_lm2$residuals) %*% GRM_oall %*% model_lm2$residuals) / sum((GRM_oall)^2)
  hat_sigma_S2 <- (sum(model_lm2$residuals^2) - hat_tau_S2 * sum(diag(GRM))) / n

  # Cross-trait variance (tau_TS and sigma_TS)
  hat_tau_TS <- as.numeric(t(model_lm2$residuals[obs_protein_index]) %*% GRM_o %*% model_lm1$residuals) / sum((GRM_o)^2)
  hat_sigma_TS <- (sum(model_lm2$residuals[obs_protein_index] * model_lm1$residuals) -
                     hat_tau_TS * sum(diag(GRM_obs_obs))) / n_obs

  # Construct Joint Covariance Matrices
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)

  I_matrix <- sparseMatrix(i = 1:min(n_obs, n), j = 1:min(n_obs, n),
                           x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n))
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix

  # Matrix Inversion (Block-wise)
  inv_Sigma22 <- matrix_inv_block(wait_matrix = Sigma22)

  # Pre-calculate conditional variance components (Schur Complement)
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- matrix_inv_Amatrix(Amatrix = A)

  # Prepare phenotype and covariate matrices for Step 2
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  X_all <- as.matrix(cbind(intercept = 1, test_df[, grepl("^X\\.", names(test_df))]))
  X_obs <- X_all[obs_protein_index, ]

  # Pre-calculate algebraic terms to optimize single-SNP testing
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))

  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))

  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all

  # Export pre-calculated parameter list
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att",
              "A22", "Atb", "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
              "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12")

  return(mget(params, envir = environment()))
}

#' Step 1 Parameter Estimation for Oracle Mixed Model (Full Observation)
#'
#' This function performs the Step 1 estimation for the Oracle model, where the
#' target phenotype is assumed to be fully observed for all samples. It estimates
#' variance components (tau and sigma) for the full cohort using a Method-of-Moments
#' approach and pre-calculates the inverse covariance matrix, residuals, and
#' projection terms for the score test.
#'
#' @param mydf A list containing the simulated data (usually the output of \code{DGP}):
#'   \itemize{
#'     \item{\code{X_all}}: Covariates for all samples.
#'     \item{\code{Y}}: The fully observed target phenotype (Oracle Y).
#'     \item{\code{GRM}}: The full Genetic Relatedness Matrix for all individuals.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) required for
#'   the Oracle score test:
#' \itemize{
#'   \item{\code{invS11_res}}: The product of the inverse covariance matrix and the
#'     null model residuals (\eqn{\Sigma^{-1} \hat{\epsilon}}).
#'   \item{\code{SX}}: The product of the inverse covariance matrix and the
#'     covariate matrix (\eqn{\Sigma^{-1} X}).
#'   \item{\code{inv_Sigma11}}: The inverse of the full covariance matrix \eqn{\Sigma}.
#'   \item{\code{XtSX_inv}}: The inverse of the information matrix for covariates
#'     (\eqn{(X^T \Sigma^{-1} X)^{-1}}).
#'   \item{\code{X}}: The covariate matrix including the intercept.
#' }
#' @export
#' @import Matrix
OracleG_typeIerror_step1 <- function(mydf) {

  # Organize full cohort data
  GRM <- mydf$GRM
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y)
  n <- nrow(test_df)

  # Prepare data for initial null model fitting
  reml_data_oracle <- data.frame(
    y = test_df$Y_all,
    x = test_df$X
  )

  # Fit a standard linear model to obtain initial residuals
  model_lm <- lm(y ~ x, data = reml_data_oracle)

  # Set diagonal of GRM to zero for Method-of-Moments variance estimation
  GRM_oall <- GRM
  diag(GRM_oall) <- 0

  # Estimate variance components (tau^2 and sigma^2) using the trace method
  # hat_tau_T2: Random effect variance component
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2)

  # hat_sigma_T2: Residual (idiosyncratic) variance component
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM))) / n

  # Construct the full cohort covariance matrix Sigma
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)

  # Compute the inverse covariance matrix using block-wise optimization
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## Pre-calculate components for the score test statistic
  # Extract outcomes and construct covariate matrix (including intercept)
  Y <- reml_data_oracle$y
  X <- cbind(intercept = 1, x = reml_data_oracle$x)

  # SX: Weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Covariance of the fixed effect estimates
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate residuals and project them using the inverse covariance matrix
  # This represents the "de-correlated" residuals used in the score numerator
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package parameters for single-SNP testing (Step 2)
  params <- c("invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  OracleG_typeIerror_step1_pars <- mget(params, envir = environment())

  return(OracleG_typeIerror_step1_pars)
}

#' Step 1 Parameter Estimation for Oracle Mixed Model (Ablation Study)
#'
#' This function performs the Step 1 estimation for the Oracle model within an
#' ablation study framework. It dynamically handles multiple covariates and
#' assumes the target phenotype is fully observed for all individuals. It
#' pre-calculates the weight matrices and residuals necessary for the score test
#' while ensuring numerical stability.
#'
#' @param mydf A list containing the simulated or real dataset:
#'   \itemize{
#'     \item{\code{X_all}}: A data frame or matrix of covariates (columns prefixed with "X.").
#'     \item{\code{Y}}: The fully observed target phenotype (Oracle Y).
#'     \item{\code{GRM}}: The full Genetic Relatedness Matrix for all individuals.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the Oracle
#' score test, optimized for multi-covariate models:
#' \itemize{
#'   \item{\code{invS11_res}}: The product of the inverse covariance matrix and
#'     the null model residuals (\eqn{\Sigma^{-1} \hat{\epsilon}}).
#'   \item{\code{SX}}: Weighted covariate matrix (\eqn{\Sigma^{-1} X}).
#'   \item{\code{inv_Sigma11}}: Inverse of the full covariance matrix \eqn{\Sigma}.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for covariates.
#'   \item{\code{X}}: The design matrix including the intercept and detected covariates.
#' }
#' @export
#' @import Matrix
OracleG_ablation_estimate <- function(mydf) {

  # Organize the full Genetic Relatedness Matrix
  GRM <- mydf$GRM

  # Prepare a data frame containing all necessary variables
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  n <- nrow(test_df)

  # Subset covariates dynamically based on the "X." prefix for the ablation study
  reml_data_oracle <- data.frame(
    y = test_df$Y_all,
    x = test_df[, grepl("^X\\.", names(test_df))]
  )

  # Construct the formula dynamically for the initial null model
  preds <- grep("^x\\.X\\.", names(reml_data_oracle), value = TRUE)
  fml <- reformulate(preds, response = "y")

  # Fit the linear null model to obtain residuals
  model_lm <- lm(fml, data = reml_data_oracle)

  # Estimate variance components (tau^2 and sigma^2) using Method-of-Moments
  # Set diagonal to zero to utilize the off-diagonal elements for tau calculation
  GRM_oall <- GRM
  diag(GRM_oall) <- 0

  # hat_tau_T2: Random effect variance component
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2)

  # hat_sigma_T2: Residual variance component with a 0.01 floor for numerical stability
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- max((a2 - hat_tau_T2 * sum(diag(GRM))) / n, 0.01)

  # Construct the full cohort covariance matrix Sigma
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)

  # Compute the inverse covariance matrix using block-wise optimization
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## Pre-calculate score test components
  # Ensure X includes an intercept for the fixed effect projection
  Y <- reml_data_oracle$y
  X <- as.matrix(cbind(intercept = 1, test_df[, grepl("^X\\.", names(test_df))]))

  # SX: The precision-weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Information matrix for fixed effects
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate de-correlated residuals (Sigma^-1 * epsilon_hat)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Return all pre-calculated terms for Step 2 GWAS testing
  params <- c("invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  OracleG_ablation_est_pars <- mget(params, envir = environment())

  return(OracleG_ablation_est_pars)
}


#' Step 1 Parameter Estimation for Observed-only Mixed Model
#'
#' This function performs the Step 1 estimation for the "Observed-only" model,
#' which restricted the analysis to individuals with non-missing target phenotypes.
#' It estimates variance components using a Method-of-Moments approach based on
#' the observed sub-GRM and pre-calculates the necessary weight matrices and
#' residuals for the subsequent score test.
#'
#' @param mydf A list containing the simulated data (usually the output of \code{DGP}):
#'   \itemize{
#'     \item{\code{X_all}}: Covariates for all samples.
#'     \item{\code{Y_obs}}: Partially observed target phenotypes (contains \code{NA}s).
#'     \item{\code{GRM}}: The full Genetic Relatedness Matrix for all individuals.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the
#' observed-only score test:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed phenotypes.
#'   \item{\code{invS11_res}}: The product of the inverse observed-covariance matrix
#'     and the null model residuals (\eqn{\Sigma_{obs}^{-1} \hat{\epsilon}_{obs}}).
#'   \item{\code{SX}}: Weighted observed covariate matrix (\eqn{\Sigma_{obs}^{-1} X_{obs}}).
#'   \item{\code{inv_Sigma11}}: Inverse of the observed-sample covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for observed covariates.
#'   \item{\code{X}}: The observed covariate design matrix including the intercept.
#' }
#' @export
#' @import Matrix
ObsG_typeIerror_step1 <- function(mydf) {

  # Identify the full Genetic Relatedness Matrix
  GRM <- mydf$GRM

  # Extract observed target phenotypes and covariates
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)

  # Obtain the indices and sub-matrix for individuals with observed data
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]

  # Prepare a clean data frame for null model fitting
  data_obs <- na.omit(data.frame(
    y = test_df$Y_obs,
    x = test_df$X
  ))
  n_obs <- nrow(data_obs)

  # Fit the initial linear null model on observed samples
  model_lm <- lm(y ~ x, data = data_obs)

  # Use Method-of-Moments for variance component estimation
  # Set diagonal to zero to utilize off-diagonal relatedness for tau estimation
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0

  # hat_tau_T2: Random effect variance component for the observed set
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2)

  # hat_sigma_T2: Residual variance component for the observed set
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs))) / n_obs

  # Construct the covariance matrix for the observed samples
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)

  # Compute the inverse matrix using block-wise optimization
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## Pre-calculate terms for the score test statistic
  # Construct the observed design matrix (intercept + covariate)
  Y <- data_obs$y
  X <- cbind(intercept = 1, x = data_obs$x)

  # SX: Precision-weighted observed covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Covariance matrix of the fixed effect estimates (observed only)
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate de-correlated residuals for the score numerator
  # This corresponds to: Sigma^-1 * (Y - X * beta_hat)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package results for single-SNP testing in Step 2
  params <- c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  ObsG_typeIerror_step1_pars <- mget(params, envir = environment())

  return(ObsG_typeIerror_step1_pars)
}



#' Step 1 Parameter Estimation for Observed-only Mixed Model
#'
#' This function performs the Step 1 estimation for the "Observed-only" model,
#' which restricts the analysis to individuals with non-missing target phenotypes.
#' It estimates variance components using a Method-of-Moments approach based on
#' the observed sub-GRM and pre-calculates the necessary weight matrices and
#' residuals for the subsequent score test.
#'
#' @param mydf A list containing the simulated data (usually the output of \code{DGP}):
#'   \itemize{
#'     \item{\code{X_all}}: Covariates for all samples.
#'     \item{\code{Y_obs}}: Partially observed target phenotypes (contains \code{NA}s).
#'     \item{\code{GRM}}: The full Genetic Relatedness Matrix for all individuals.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the
#' observed-only score test:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed phenotypes.
#'   \item{\code{invS11_res}}: The product of the inverse observed-covariance matrix
#'     and the null model residuals (\eqn{\Sigma_{obs}^{-1} \hat{\epsilon}_{obs}}).
#'   \item{\code{SX}}: Weighted observed covariate matrix (\eqn{\Sigma_{obs}^{-1} X_{obs}}).
#'   \item{\code{inv_Sigma11}}: Inverse of the observed-sample covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for observed covariates.
#'   \item{\code{X}}: The observed covariate design matrix including the intercept.
#' }
#' @export
#' @import Matrix
ObsG_typeIerror_step1 <- function(mydf) {

  # Identify the full Genetic Relatedness Matrix
  GRM <- mydf$GRM

  # Extract observed target phenotypes and covariates
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)

  # Obtain the indices and sub-matrix for individuals with observed data
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]

  # Prepare a clean data frame for null model fitting
  # We use only the subset of samples where the protein is measured
  data_obs <- na.omit(data.frame(
    y = test_df$Y_obs,
    x = test_df$X
  ))
  n_obs <- nrow(data_obs)

  # Fit the initial linear null model on observed samples to get OLS residuals
  model_lm <- lm(y ~ x, data = data_obs)

  # Method-of-Moments (MoM) for variance component estimation
  # Sigma = tau^2 * GRM + sigma^2 * I
  # We use the off-diagonal trace of residuals to solve for tau^2
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0

  # Quadratic form: e' * (GRM - diag(GRM)) * e
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2)

  # Solve for the idiosyncratic error variance sigma^2
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs))) / n_obs

  # Construct the estimated V-matrix (Covariance Matrix)
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)

  # Efficient block-wise inversion
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## Pre-calculate terms for the score test statistic numerator and variance
  # Design matrix includes intercept for correct fixed-effect projection
  Y <- data_obs$y
  X <- cbind(intercept = 1, x = data_obs$x)

  # Projecting covariates into the precision space
  SX <- inv_Sigma11 %*% X

  # Information matrix inverse: (X' Sigma^-1 X)^-1
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate decorrelated residuals: Sigma^-1 * (Y - X * beta_GLS)
  # This is the "Score" part of the numerator
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package results for Step 2 GWAS loop
  params <- c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  ObsG_typeIerror_step1_pars <- mget(params, envir = environment())

  return(ObsG_typeIerror_step1_pars)
}


#' Step 1 Parameter Estimation for Observed-only Mixed Model (Ablation Study)
#'
#' This function performs the Step 1 estimation for the Observed-only model
#' within an ablation study framework. It dynamically handles multiple
#' covariates and restricts the analysis to individuals with non-missing target
#' phenotypes. It estimates variance components via Method-of-Moments and
#' pre-calculates score test components while ensuring numerical stability.
#'
#' @param mydf A list containing the simulated or real dataset:
#'   \itemize{
#'     \item{\code{X_all}}: A data frame or matrix of covariates (columns prefixed with "X.").
#'     \item{\code{Y_obs}}: Partially observed target phenotypes.
#'     \item{\code{GRM}}: The full Genetic Relatedness Matrix.
#'   }
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) optimized for
#' ablation study multi-covariate models:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed phenotypes.
#'   \item{\code{invS11_res}}: The product of the inverse observed-covariance matrix
#'     and the null model residuals (\eqn{\Sigma_{obs}^{-1} \hat{\epsilon}_{obs}}).
#'   \item{\code{SX}}: Weighted observed covariate matrix (\eqn{\Sigma_{obs}^{-1} X_{obs}}).
#'   \item{\code{inv_Sigma11}}: Inverse of the observed-sample covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for observed covariates.
#'   \item{\code{X}}: The observed design matrix including the intercept and detected covariates.
#' }
#' @export
#' @import Matrix
ObsG_ablation_estimate <- function(mydf) {

  # Identify the Genetic Relatedness Matrix
  GRM <- mydf$GRM

  # Extract data components
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)

  # Subset indices and GRM for individuals with observed phenotypes
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]

  # Prepare clean data frame for ablation study null model
  # Dynamically capture covariates prefixed with "X."
  data_obs <- na.omit(data.frame(
    y = test_df$Y_obs,
    test_df[, grepl("^X\\.", names(test_df))]
  ))
  n_obs <- nrow(data_obs)

  # Dynamically construct the formula for the null model fit
  preds <- grep("^X\\.", names(data_obs), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm <- lm(fml, data = data_obs)

  # Method-of-Moments Variance Component Estimation
  # Use off-diagonal elements for random effect variance (tau^2)
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0

  # hat_tau_T2: Random effect variance component
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2)

  # hat_sigma_T2: Residual variance with 0.01 floor protection for numerical stability
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- max((a2 - hat_tau_T2 * sum(diag(GRM_obs_obs))) / n_obs, 0.01)

  # Construct and invert the observed-sample covariance matrix
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## Pre-calculate terms for the Score Test statistic
  # Ensure X includes an intercept and all relevant covariates from the data_obs subset
  Y <- data_obs$y
  # Note: grepl matches the "X." prefix in the columns of the subsetted data_obs
  X <- as.matrix(cbind(intercept = 1, data_obs[, grepl("^X\\.", names(data_obs))]))

  # SX: Precision-weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Information matrix inverse for observed fixed effects
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate decorrelated residuals: Sigma^-1 * (Y - X * beta_GLS)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package and return results
  params <- c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  ObsG_ablation_est_pars <- mget(params, envir = environment())

  return(ObsG_ablation_est_pars)
}


#' Step 1 Parameter Estimation for Syn-PALM (Independent Samples)
#'
#' This function performs Step 1 estimation for the Syn-PALM method assuming
#' independent individuals (no relatedness). It estimates variance and covariance
#' components between the target protein and surrogate phenotype using standard
#' OLS residuals, and pre-calculates the matrices required for the score test.
#' It is specifically designed for use with a subset of independent samples.
#'
#' @param mydf A list containing the full simulated data (X_all, Y, Y_obs, S).
#' @param independent_indices A numeric vector of indices specifying which
#'   individuals to include (the "independent" subset).
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) optimized for
#' independent samples:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed
#'     phenotypes within the independent subset.
#'   \item{\code{V11}}: Inverse of the conditional variance matrix (Schur Complement).
#'   \item{\code{B_mat}}: Projection matrix for beta estimation under GLS.
#'   \item{\code{Sigma12_Sigma22inv}}: Projection of surrogates onto the target space.
#'   \item{... and other matrices required for Step 2 GWAS scanning.}
#' }
#' @export
#' @import Matrix
SynSurr_typeIerror_step1 <- function(mydf, independent_indices) {

  # Subset data for independent samples
  test_df_independent <- data.frame(
    X = mydf$X_all[independent_indices],
    Y_all = mydf$Y[independent_indices],
    Y_obs = mydf$Y_obs[independent_indices],
    yhat = mydf$S[independent_indices]
  )

  n <- nrow(test_df_independent)

  # Identify observed indices within the subset
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  n_obs <- length(obs_protein_index)

  # Prepare datasets for variance estimation
  reml_data_SynSurr <- na.omit(data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  ))

  lm_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  )

  # Fit null models via OLS to get residuals
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurr)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurr)

  # Estimate variance/covariance components (Scalar values for independent samples)
  hat_sigma_T2 <- var(model_lm1$residuals)
  hat_sigma_S2 <- var(model_lm2$residuals)
  hat_sigma_TS <- cov(model_lm2$residuals[obs_protein_index], model_lm1$residuals)

  # Construct Diagonal Covariance Matrices
  # Sigma11: Observed Target Variance
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)

  # Sigma12: Cross-trait Covariance (Target vs Surrogate)
  # Represents the identity-like mapping of errors between traits
  I_matrix <- sparseMatrix(
    i = 1:min(n_obs, n), j = 1:min(n_obs, n),
    x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n)
  )
  Sigma12 <- I_matrix

  # Sigma22: Full Surrogate Variance
  Sigma22 <- hat_sigma_S2 * Diagonal(n)

  # Matrix Inversion via Cholesky Decomposition (more efficient for independent data)
  inv_Sigma11 <- chol2inv(chol(Sigma11))
  inv_Sigma22 <- chol2inv(chol(Sigma22))

  # Pre-calculate SynSurr components (Schur Complement logic)
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- chol2inv(chol(A))

  ## Pre-calculate GLS estimation components
  Y <- test_df_independent$Y_obs[obs_protein_index]
  hatY <- test_df_independent$yhat

  # Design matrices (including intercept)
  X_all <- as.matrix(cbind(intercept = 1, x = test_df_independent$X))
  X_obs <- X_all[obs_protein_index, ]

  # Algebraic terms for Step 2 GWAS loop speed-up
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))

  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))

  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all

  # Package results
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att",
              "A22", "Atb", "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
              "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12")

  return(mget(params, envir = environment()))
}

#' Step 1 Parameter Estimation for Oracle Model (Independent Samples)
#'
#' This function performs Step 1 estimation for the Oracle model assuming
#' independent individuals. It calculates the variance of the target phenotype
#' using the full cohort (all samples observed) and pre-calculates the weight
#' matrices and residuals required for the score test. This serves as the
#' theoretical gold standard for efficiency and power in independent sample simulations.
#'
#' @param mydf A list containing the full simulated data (X_all, Y).
#' @param independent_indices A numeric vector of indices specifying the
#'   independent subset of individuals to be analyzed.
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the Oracle
#' score test:
#' \itemize{
#'   \item{\code{invS11_res}}: The product of the inverse covariance matrix and
#'     the full-sample null residuals (\eqn{\Sigma^{-1} \hat{\epsilon}}).
#'   \item{\code{SX}}: Weighted covariate matrix (\eqn{\Sigma^{-1} X}).
#'   \item{\code{inv_Sigma11}}: Inverse of the diagonal covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for fixed effects.
#'   \item{\code{X}}: The covariate design matrix including the intercept.
#' }
#' @export
#' @import Matrix
Oracle_typeIerror_step1 <- function(mydf, independent_indices) {

  # Subset the data for the specified independent individuals
  test_df_independent <- data.frame(
    X = mydf$X_all[independent_indices],
    Y_all = mydf$Y[independent_indices],
    Y_obs = mydf$Y_obs[independent_indices],
    yhat = mydf$S[independent_indices]
  )

  n <- nrow(test_df_independent)

  # Prepare data for the Oracle null model (assuming Y is fully observed)
  data_oracle <- data.frame(
    y = test_df_independent$Y_all,
    x = test_df_independent$X
  )

  # Fit a standard linear model to estimate residual variance
  model_lm <- lm(y ~ x, data = data_oracle)

  # Estimate the error variance (scalar, as samples are independent)
  hat_sigma_T2 <- var(model_lm$residuals)

  # Construct the diagonal covariance matrix Sigma
  Sigma11 <- hat_sigma_T2 * Diagonal(n)

  # Invert Sigma using Cholesky decomposition for efficiency
  chol_Sigma11 <- chol(Sigma11)
  inv_Sigma11 <- chol2inv(chol_Sigma11)

  ## Pre-calculate terms for the single-SNP score test (Step 2)
  # Ensure X includes an intercept for the fixed effect projection
  Y <- data_oracle$y
  X <- cbind(intercept = 1, x = data_oracle$x)

  # SX: Precision-weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Variance-covariance matrix of the fixed effect estimates
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate decorrelated residuals: Sigma^-1 * (Y - X * beta_hat)
  # This serves as the adjusted outcome vector in the score test numerator
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package and return parameters
  params <- c("invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  Oracle_typeIerror_step1_pars <- mget(params, envir = environment())

  return(Oracle_typeIerror_step1_pars)
}

#' Step 1 Parameter Estimation for Observed-only Model (Independent Samples)
#'
#' This function performs Step 1 estimation for the "Observed-only" model
#' assuming independent individuals. It restricts the analysis to the subset of
#' individuals with non-missing target phenotypes within a specified independent
#' group. It estimates the residual variance via OLS and pre-calculates the
#' weight matrices and residuals needed for the score test.
#'
#' @param mydf A list containing the full simulated data (X_all, Y_obs).
#' @param independent_indices A numeric vector of indices specifying the
#'   independent subset of individuals to be analyzed.
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the
#' observed-only score test on independent samples:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed
#'     phenotypes relative to the independent subset.
#'   \item{\code{invS11_res}}: The product of the inverse observed-covariance matrix
#'     and the null residuals (\eqn{\Sigma_{obs}^{-1} \hat{\epsilon}_{obs}}).
#'   \item{\code{SX}}: Weighted observed covariate matrix (\eqn{\Sigma_{obs}^{-1} X_{obs}}).
#'   \item{\code{inv_Sigma11}}: Inverse of the diagonal observed covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for observed fixed effects.
#'   \item{\code{X}}: The observed covariate design matrix including the intercept.
#' }
#' @export
#' @import Matrix
Obs_typeIerror_step1 <- function(mydf, independent_indices) {

  # Subset the data for the specified independent individuals
  test_df_independent <- data.frame(
    X = mydf$X_all[independent_indices],
    Y_all = mydf$Y[independent_indices],
    Y_obs = mydf$Y_obs[independent_indices],
    yhat = mydf$S[independent_indices]
  )

  # Identify indices of individuals with observed phenotypes within the subset
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))

  # Prepare clean data frame for null model fitting (observed only)
  data_obs <- na.omit(data.frame(
    y = test_df_independent$Y_obs,
    x = test_df_independent$X
  ))
  n_obs <- nrow(data_obs)

  # Fit a standard linear model to estimate residual variance for observed samples
  model_lm <- lm(y ~ x, data = data_obs)

  # Estimate error variance (scalar, as samples are independent)
  hat_sigma_T2 <- var(model_lm$residuals)

  # Construct the diagonal covariance matrix for observed samples
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)

  # Efficiently invert the diagonal matrix using Cholesky
  chol_Sigma11 <- chol(Sigma11)
  inv_Sigma11 <- chol2inv(chol_Sigma11)

  ## Pre-calculate score test components
  # Construct the design matrix (intercept + covariate) for the observed subset
  Y <- data_obs$y
  X <- cbind(intercept = 1, x = data_obs$x)

  # SX: Precision-weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Information matrix inverse for observed fixed effects
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate decorrelated residuals: Sigma^-1 * (Y - X * beta_hat)
  # Represents the score part of the statistic numerator
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package results for Step 2 GWAS loop
  params <- c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  Obs_typeIerror_step1_pars <- mget(params, envir = environment())

  return(Obs_typeIerror_step1_pars)
}

#' Step 1 Parameter Estimation for Syn-PALM Ablation Study (Independent Samples)
#'
#' This function performs Step 1 estimation for ablation studies within the
#' Syn-PALM framework, specifically tailored for independent individuals. It
#' dynamically detects multiple covariates and assumes a diagonal covariance
#' structure (no relatedness). It pre-calculates the necessary conditional
#' variance and projection matrices using high-speed Cholesky decomposition.
#'
#' @param mydf A list containing the full dataset:
#'   \itemize{
#'     \item{\code{X_all}}: A data frame or matrix of covariates (columns prefixed with "X.").
#'     \item{\code{Y}}: Oracle phenotypes.
#'     \item{\code{Y_obs}}: Partially observed target phenotypes.
#'     \item{\code{S}}: Synthetic/Surrogate phenotypes.
#'   }
#' @param independent_indices A numeric vector of indices specifying the
#'   independent subset of individuals for the ablation analysis.
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) optimized for
#'   multi-covariate models without relatedness:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed
#'     phenotypes within the independent subset.
#'   \item{\code{V11}}: Inverse of the conditional variance matrix (Schur Complement).
#'   \item{\code{B_mat}}: Projection matrix for fixed-effect estimation.
#'   \item{... and other algebraic terms for efficient single-SNP testing.}
#' }
#' @export
#' @import Matrix
SynSurr_ablation_estimate <- function(mydf, independent_indices) {

  # Subset data for independent samples
  test_df_independent <- data.frame(
    X = mydf$X_all[independent_indices, ],
    Y_all = mydf$Y[independent_indices],
    Y_obs = mydf$Y_obs[independent_indices],
    yhat = mydf$S[independent_indices]
  )

  n <- nrow(test_df_independent)

  # Identify observed indices within the subset
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  n_obs <- length(obs_protein_index)

  # Prepare datasets for variance estimation (dynamic covariate selection)
  reml_data_SynSurr <- na.omit(data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    test_df_independent[, grepl("^X\\.", names(test_df_independent))]
  ))

  lm_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    test_df_independent[, grepl("^X\\.", names(test_df_independent))]
  )

  # Construct formulas dynamically for Y and Surrogate S
  preds_y <- grep("^x\\.X\\.", names(reml_data_SynSurr), value = TRUE)
  fml_y <- reformulate(preds_y, response = "y")
  model_lm1 <- lm(fml_y, data = reml_data_SynSurr)

  preds_s <- grep("^x\\.X\\.", names(lm_data_SynSurr), value = TRUE)
  fml_s <- reformulate(preds_s, response = "haty")
  model_lm2 <- lm(fml_s, data = lm_data_SynSurr)

  # Variance/Covariance components (Scalar estimates for independent samples)
  hat_sigma_T2 <- var(model_lm1$residuals)
  hat_sigma_S2 <- var(model_lm2$residuals)
  hat_sigma_TS <- cov(model_lm2$residuals[obs_protein_index], model_lm1$residuals)

  # Construct Diagonal Covariance Matrices
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  Sigma22 <- hat_sigma_S2 * Diagonal(n)

  # Sigma12: Cross-trait Covariance mapping
  I_matrix <- sparseMatrix(
    i = 1:min(n_obs, n), j = 1:min(n_obs, n),
    x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n)
  )
  Sigma12 <- I_matrix

  # Matrix Inversion via Cholesky (optimized for Diagonal matrices)
  inv_Sigma11 <- chol2inv(chol(Sigma11))
  inv_Sigma22 <- chol2inv(chol(Sigma22))

  # Schur Complement and precision matrix V11 calculation
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- chol2inv(chol(A))

  # Pre-calculate GLS and projection components
  Y <- test_df_independent$Y_obs[obs_protein_index]
  hatY <- test_df_independent$yhat

  # Design matrices with intercept and dynamic covariates
  X_all <- as.matrix(cbind(intercept = 1,
                           test_df_independent[, grepl("^X\\.", names(test_df_independent))]))
  X_obs <- X_all[obs_protein_index, ]

  # Speed up Step 2 by pre-calculating constant matrix products
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))

  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))

  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all

  # Export parameter set
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att",
              "A22", "Atb", "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
              "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12")

  return(mget(params, envir = environment()))
}

#' Step 1 Parameter Estimation for Observed-only Ablation Study (Independent Samples)
#'
#' This function performs Step 1 estimation for the "Observed-only" model within
#' an ablation study framework for independent individuals. It dynamically
#' handles multiple covariates (prefixed with "X.") and restricts analysis to
#' samples with non-missing target phenotypes. It serves as a baseline to
#' evaluate the power gain of Syn-PALM under different covariate settings.
#'
#' @param mydf A list containing the full dataset:
#'   \itemize{
#'     \item{\code{X_all}}: A data frame or matrix of covariates (columns prefixed with "X.").
#'     \item{\code{Y_obs}}: Partially observed target phenotypes.
#'   }
#' @param independent_indices A numeric vector of indices specifying the
#'   independent subset of individuals for the ablation analysis.
#'
#' @return A list of pre-calculated parameters (\code{step1_pars}) for the
#' observed-only score test, optimized for independent samples:
#' \itemize{
#'   \item{\code{obs_protein_index}}: Indices of individuals with observed
#'     phenotypes relative to the independent subset.
#'   \item{\code{invS11_res}}: The product of the inverse observed-covariance matrix
#'     and the null residuals (\eqn{\Sigma_{obs}^{-1} \hat{\epsilon}_{obs}}).
#'   \item{\code{SX}}: Weighted observed covariate matrix (\eqn{\Sigma_{obs}^{-1} X_{obs}}).
#'   \item{\code{inv_Sigma11}}: Inverse of the diagonal observed covariance matrix.
#'   \item{\code{XtSX_inv}}: Inverse of the information matrix for observed fixed effects.
#'   \item{\code{X}}: The observed design matrix including the intercept and detected covariates.
#' }
#' @export
#' @import Matrix
Obs_ablation_estimate <- function(mydf, independent_indices) {

  # Subset data for the specified independent individuals
  test_df_independent <- data.frame(
    X = mydf$X_all[independent_indices, ],
    Y_all = mydf$Y[independent_indices],
    Y_obs = mydf$Y_obs[independent_indices],
    yhat = mydf$S[independent_indices]
  )

  # Identify indices of individuals with observed phenotypes
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))

  # Prepare clean data frame for ablation study null model (observed only)
  # Dynamically capture covariates prefixed with "X."
  data_obs <- na.omit(data.frame(
    y = test_df_independent$Y_obs,
    test_df_independent[, grepl("^X\\.", names(test_df_independent))]
  ))
  n_obs <- nrow(data_obs)

  # Dynamically construct the formula for the null model fit
  preds <- grep("^X\\.", names(data_obs), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm <- lm(fml, data = data_obs)

  # Estimate residual variance (scalar for independent samples)
  hat_sigma_T2 <- var(model_lm$residuals)

  # Construct diagonal covariance matrix for observed samples
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)

  # Efficiently invert the diagonal matrix via Cholesky
  chol_Sigma11 <- chol(Sigma11)
  inv_Sigma11 <- chol2inv(chol_Sigma11)

  ## Pre-calculate terms for Step 2 GWAS Score Test
  # Design matrix includes intercept and all dynamic covariates
  Y <- data_obs$y
  X <- as.matrix(cbind(intercept = 1, data_obs[, grepl("^X\\.", names(data_obs))]))

  # SX: Precision-weighted covariate matrix
  SX <- inv_Sigma11 %*% X

  # XtSX_inv: Information matrix inverse for observed fixed effects
  XtSX_inv <- solve(crossprod(X, SX))

  # Calculate decorrelated residuals: Sigma^-1 * (Y - X * beta_hat)
  # This corresponds to the score part of the numerator in the GWAS loop
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  # Package and return results
  params <- c("obs_protein_index", "invS11_res", "SX", "inv_Sigma11", "XtSX_inv", "X")
  Obs_ablation_est_pars <- mget(params, envir = environment())

  return(Obs_ablation_est_pars)
}

Obs_ablation_estimate <- function(mydf,independent_indices) {

  test_df_independent <- data.frame(X = mydf$X_all[independent_indices,],
                                    Y_all = mydf$Y[independent_indices],
                                    Y_obs = mydf$Y_obs[independent_indices],
                                    yhat = mydf$S[independent_indices])

  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))

  data_obs <- data.frame(
    y = test_df_independent$Y_obs,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent))]
  )%>% na.omit()

  n_obs <- nrow(data_obs)

  preds <- grep("^x\\.X\\.", names(data_obs), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm <- lm(fml, data = data_obs)

  hat_sigma_T2 <- var(model_lm$residuals)
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)

  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11)


  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- as.matrix(cbind(rep(1,length(Y)),data_obs[, grepl("^x\\.X\\.", names(data_obs))]))

  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual

  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X"
  )

  Obs_ablation_est_pars <- mget(params, envir = environment())

  return(Obs_ablation_est_pars)

}


#' Integrated Syn-PALM Score Test for Power Simulation (Single SNP)
#'
#' This high-level wrapper function performs the complete Syn-PALM workflow
#' for a single genetic variant. It estimates variance components using
#' computationally efficient methods, constructs the necessary block-covariance
#' matrices (Sigma11, Sigma12, Sigma22), and executes the score test. It is
#' specifically optimized for large-scale power simulations where REML
#' convergence is too slow.
#'
#' @param mydf A list object containing the simulated data:
#'   \itemize{
#'     \item{\code{Y_obs}}: Partially observed target phenotype.
#'     \item{\code{S}}: Synthetic/Surrogate phenotype.
#'     \item{\code{X_all}}: Covariates.
#'     \item{\code{G_all}}: Genotype vector of the candidate SNP.
#'     \item{\code{GRM}}: Genetic Relatedness Matrix.
#'   }
#'
#' @return A list containing the results of the score test:
#' \itemize{
#'   \item{\code{score}}: The calculated score statistic.
#'   \item{\code{p_value}}: The p-value for the association between G and Y.
#'   \item{... and other statistics from \code{score_test_SynSurrG_single}.}
#' }
#'
#' @details
#' The function operates in three phases:
#' \enumerate{
#'   \item \strong{Variance Estimation}: Uses a non-iterative approach (similar to
#'         H-E regression) to estimate trait-specific and cross-trait variance
#'         components (\eqn{\tau^2_T, \sigma^2_T, \tau^2_S, \sigma^2_S, \tau_{TS}, \sigma_{TS}}).
#'   \item \strong{Step 1 (Pre-calculation)}: Constructs the joint covariance structure
#'         and performs block-wise matrix inversion and Schur complement calculation
#'         to derive the precision matrix \eqn{V_{11}}.
#'   \item \strong{Step 2 (Testing)}: Calls the single-SNP score test function
#'         using the pre-calculated parameters.
#' }
#' @export
#' @import Matrix
score_test_SynSurrG_power_sim_one_g_to_one_Y <- function(mydf) {

  # Setup local data frame for modeling
  test_df <- data.frame(
    X = mydf$X_all,
    Y_all = mydf$Y,
    Y_obs = mydf$Y_obs,
    yhat = mydf$S
  )

  n <- nrow(test_df)
  GRM <- mydf$GRM

  # Identify observed indices and subset the GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index, 1:n]
  n_obs <- length(obs_protein_index)

  # Initial OLS to obtain residuals for variance component estimation
  reml_data_SynSurrG <- na.omit(data.frame(y = test_df$Y_obs, haty = test_df$yhat, x = test_df$X))
  lm_data_SynSurrG <- data.frame(y = test_df$Y_obs, haty = test_df$yhat, x = test_df$X)

  model_lm1 <- lm(y ~ x, data = reml_data_SynSurrG)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurrG)

  # Variance Component Estimation (Fast H-E style approach)
  # Removes diagonal for off-diagonal relatedness calculations
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0

  # Estimate tau^2 and sigma^2 for Target Trait (T)
  a1 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2)
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs))) / nrow(GRM_obs_obs)

  # Estimate tau^2 and sigma^2 for Surrogate Trait (S)
  a1 <- as.numeric(t(model_lm2$residuals) %*% GRM_oall %*% model_lm2$residuals)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2)
  a2 <- sum(model_lm2$residuals * model_lm2$residuals)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM))) / nrow(GRM)

  # Estimate Cross-trait variance components (TS)
  a1 <- as.numeric(t(model_lm2$residuals[obs_protein_index]) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_TS <- a1 / sum((GRM_o)^2)
  a2 <- sum(model_lm2$residuals[obs_protein_index] * model_lm1$residuals)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs))) / nrow(GRM_obs_obs)

  # Construct the Joint Covariance Blocks
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)

  I_matrix <- sparseMatrix(
    i = 1:min(n_obs, n), j = 1:min(n_obs, n),
    x = rep(hat_sigma_TS, min(n_obs, n)), dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)

  # Inversions and Schur Complement (using block-wise optimization)
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)
  inv_Sigma22 <- matrix_inv_block(wait_matrix = Sigma22)

  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  V11 <- matrix_inv_Amatrix(Amatrix = A)

  # Fixed Effects Pre-calculation (GLS approach)
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  X_all <- as.matrix(cbind(intercept = 1, test_df$X))
  X_obs <- X_all[obs_protein_index, ]

  # Prepare intermediate matrices for the score test
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1)
  B_mat <- solve(Btt2, t(Btt1))
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  bt2 <- Sigma12_Sigma22inv %*% X_all

  # Consolidate parameters for Step 2
  params <- c("obs_protein_index", "V11", "Y", "hatY", "X_obs", "Att",
              "A22", "Atb", "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
              "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12")
  SynSurrG_step1_pars <- mget(params, envir = environment())

  # Execute Step 2 Score Test
  G_all <- mydf$G_all
  results <- score_test_SynSurrG_single(G_all, step1_pars = SynSurrG_step1_pars)

  return(results)
}

#' Integrated Observed-only Score Test for Power Simulation (Baseline Mixed Model)
#'
#' This function performs a standard mixed-model score test using only individuals
#' with observed target phenotypes. It functions as a baseline comparator in power
#' simulations. Like its Syn-PALM counterpart, it utilizes a fast H-E regression
#' approach for variance component estimation to ensure computational feasibility
#' across thousands of simulation replicates.
#'
#' @param mydf A list object containing the simulated data:
#'   \itemize{
#'     \item{\code{Y_obs}}: Partially observed target phenotype.
#'     \item{\code{X_all}}: Covariates (e.g., PCs, Age).
#'     \item{\code{G_all}}: Genotype vector of the candidate SNP.
#'     \item{\code{GRM}}: Genetic Relatedness Matrix.
#'   }
#'
#' @return A named numeric vector:
#' \itemize{
#'   \item{\code{T_score_ObsG}}: The calculated score statistic for the observed-only model.
#'   \item{\code{negative_log10_pval_ObsG}}: The -log10 p-value of the association.
#' }
#'
#' @details
#' The workflow follows standard GWAS mixed-model principles:
#' \enumerate{
#'   \item \strong{Variance Component Estimation}: Estimates \eqn{\tau^2_T} (genetic variance)
#'         and \eqn{\sigma^2_T} (residual variance) using the H-E method on the
#'         observed subset.
#'   \item \strong{Precision Calculation}: Construct and invert the subset covariance
#'         matrix \eqn{\Sigma_{11} = \tau^2_T K_{obs} + \sigma^2_T I}.
#'   \item \strong{Score Statistic}: Calculates the score \eqn{U = G_{obs}' \Sigma_{11}^{-1} \hat{\epsilon}}
#'         and its variance \eqn{V(U)}, accounting for the uncertainty in fixed effects.
#' }
#' @export
#' @import Matrix
score_test_ObsG_power_sim_one_g_to_one_Y <- function(mydf) {

  # Local data setup
  test_df <- data.frame(
    X = mydf$X_all,
    Y_all = mydf$Y,
    Y_obs = mydf$Y_obs,
    yhat = mydf$S
  )

  n <- nrow(test_df)
  GRM <- mydf$GRM

  # Subset GRM and identify observed individuals
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  n_obs <- length(obs_protein_index)

  # Initial OLS on observed data to get residuals for H-E estimation
  data_obsG <- na.omit(data.frame(y = test_df$Y_obs, x = test_df$X))
  model_lm1 <- lm(y ~ x, data = data_obsG)

  # Fast Variance Component Estimation (H-E Method)
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0 # Mask diagonal for relatedness variance calculation

  # Estimate tau^2 (Random effect variance)
  a1 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2)

  # Estimate sigma^2 (Residual variance)
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs))) / nrow(GRM_obs_obs)

  # Construct observed-only covariance matrix
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)

  # Block-wise matrix inversion for the observed subset
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## GLS Fixed Effect Estimation and Residual Calculation
  Y <- test_df$Y_obs[obs_protein_index]
  X_all <- as.matrix(cbind(intercept = 1, test_df$X))
  X_obs <- X_all[obs_protein_index, ]
  G_all <- mydf$G_all

  # Weighted design matrix and residuals
  SX <- inv_Sigma11 %*% X_obs
  XtSX_inv <- solve(crossprod(X_obs, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X_obs %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual

  ## Score Test Statistic Construction
  G_obs <- G_all[obs_protein_index]

  # Score Numerator (SU)
  SU <- as.numeric(crossprod(G_obs, invS11_res))

  # Score Variance (VU) - Adjusted for fixed effects uncertainty
  A <- crossprod(SX, G_obs)
  VU <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * (XtSX_inv %*% A)))

  # Result Calculation (Chi-square test)
  results <- compute_score(SU, VU)
  names(results) <- c("T_score_ObsG", "negative_log10_pval_ObsG")

  return(results)
}

#' Integrated Oracle Score Test for Power Simulation (Ideal Case Mixed Model)
#'
#' This function performs a mixed-model score test assuming the target phenotype
#' (Y) is fully observed for the entire cohort. It serves as the "Gold Standard"
#' or theoretical upper bound in power simulations, allowing you to quantify
#' how much information Syn-PALM recovers relative to the ideal full-data scenario.
#'
#' @param mydf A list object containing the simulated data:
#'   \itemize{
#'     \item{\code{Y}}: The Oracle (fully observed) target phenotype.
#'     \item{\code{X_all}}: Covariates for the full cohort.
#'     \item{\code{G_all}}: Genotype vector for the full cohort.
#'     \item{\code{GRM}}: Genetic Relatedness Matrix for the full cohort.
#'   }
#'
#' @return A named numeric vector:
#' \itemize{
#'   \item{\code{T_score_OracleG}}: The calculated score statistic for the Oracle model.
#'   \item{\code{negative_log10_pval_OracleG}}: The -log10 p-value of the association.
#' }
#'
#' @details
#' The function follows the same computational logic as the observed-only and
#' Syn-PALM models to ensure comparability:
#' \enumerate{
#'   \item \strong{Variance Estimation}: Estimates \eqn{\tau^2_T} and \eqn{\sigma^2_T}
#'         using the H-E method applied to the full \eqn{N \times N} GRM.
#'   \item \strong{Precision Calculation}: Inverts the full covariance matrix
#'         \eqn{\Sigma = \tau^2_T K + \sigma^2_T I}.
#'   \item \strong{Score Test}: Computes the association between \code{G_all} and
#'         \code{Y} using the full sample size \eqn{N}.
#' }
#' @export
#' @import Matrix
score_test_OracleG_power_sim_one_g_to_one_Y <- function(mydf) {

  # Local data setup using the Oracle (complete) phenotype
  test_df <- data.frame(
    X = mydf$X_all,
    Y_all = mydf$Y,
    Y_obs = mydf$Y_obs,
    yhat = mydf$S
  )

  n <- nrow(test_df)
  GRM <- mydf$GRM
  G_all <- mydf$G_all

  # Fit null model on complete data to obtain residuals
  lm_data_OracleG <- data.frame(y = test_df$Y_all, x = test_df$X)
  model_lm <- lm(y ~ x, data = lm_data_OracleG)

  # Fast Variance Component Estimation (H-E Method for full cohort)
  GRM_oall <- GRM
  diag(GRM_oall) <- 0 # Mask diagonal

  # Estimate tau^2 (Random effect variance)
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2)

  # Estimate sigma^2 (Residual variance)
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM))) / nrow(GRM)

  # Construct full-sample covariance matrix
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)

  # Invert full covariance matrix via block-wise optimization
  inv_Sigma11 <- matrix_inv_block(wait_matrix = Sigma11)

  ## GLS Fixed Effect Estimation
  Y <- lm_data_OracleG$y
  X <- as.matrix(cbind(intercept = 1, lm_data_OracleG$x))

  # Calculate weighted design matrix and residuals
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual

  ## Oracle Score Test Statistic
  # Score Numerator (SU)
  SU <- as.numeric(crossprod(G_all, invS11_res))

  # Score Variance (VU) adjusted for fixed effects
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))

  # Final results (Chi-square test)
  results <- compute_score(SU, VU)
  names(results) <- c("T_score_OracleG", "negative_log10_pval_OracleG")

  return(results)
}


#' Reset and Recreate a Temporary Directory
#'
#' This function checks if a directory exists, removes it if it does,
#' and then creates a fresh one.
#'
#' @param dir Character string. The path to the directory to be reset.
#' @export

reset_temp_dir <- function(dir) {
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

#' Create a Sample File in BGEN/QCTOOL Format
#'
#' This function generates a \code{.sample} file with the required two-line header,
#' compatible with UKB and other genomic tools.
#'
#' @param eid Vector of sample IDs (e.g., UKB eids).
#' @param covariates A data.frame with the same number of rows as \code{eid}, containing additional covariates.
#' @param file_out Character string for the output filename.
#'
#' @return None. Writes a file to the specified path.
#' @export
make_sample_file <- function(eid, covariates = NULL, file_out = "output.sample") {
  # eid: vector (UKB ID or other sample ID)
  # covariates: data.frame, nrow = length(eid), columns are additional covariates
  # file_out: output filename

  n <- length(eid)

  # Basic columns
  df <- data.frame(
    ID_1 = eid,
    ID_2 = eid,
    missing = 0
  )

  # Add covariates
  if (!is.null(covariates)) {
    stopifnot(nrow(covariates) == n)
    df <- cbind(df, covariates)
  }

  # === Write file ===
  # First line: column names
  cat(paste(colnames(df), collapse = " "), "\n", file = file_out)

  # Second line: type indicators
  #   ID_1, ID_2, missing -> "0"
  #   Other columns -> "D" (stands for Discrete/covariate)
  type_line <- c("0", "0", "0", rep("D", ncol(df) - 3))
  cat(paste(type_line, collapse = " "), "\n", file = file_out, append = TRUE)

  # Subsequent lines: sample information
  write.table(
    df,
    file = file_out,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )

  message("✅ .sample file written to: ", file_out)
}

#' Split Samples into Training and Testing Sets by Blocks
#'
#' This function selects a subset of blocks to form a training set such that the total
#' number of samples is as close as possible to a specified proportion of the total
#' sample size. It uses a random search with a greedy approach.
#'
#' @param hh A list where each element is a numeric vector of sample indices belonging to a block.
#' @param prop Numeric. The target proportion for the training set. Default is 0.8.
#' @param iterations Integer. The number of random search iterations to find the optimal combination. Default is 5000.
#' @param seed Integer. Optional random seed for reproducibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item{train_ids}: Vector of sample IDs assigned to the training set.
#'   \item{test_ids}: Vector of sample IDs assigned to the test set.
#'   \item{split}: A character vector indexed by sample ID indicating "train" or "test".
#'   \item{sizes}: A named vector with counts for train, test, and total samples.
#'   \item{target}: A named vector showing the target size, achieved size, and the difference.
#'   \item{block_sizes}: A vector containing the size of each input block.
#' }
#' @export
split_by_blocks <- function(hh, prop = 0.8, iterations = 5000, seed = NULL) {
  # hh: list, each element is a vector of sample indices for a block (e.g., hh[[1]], hh[[2]]...)
  # prop: training set proportion (default 0.8)
  # iterations: number of random search iterations; larger values yield better approximations (1e3~1e4 is usually sufficient)
  # seed: optional, set random seed for reproducibility

  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.list(hh), length(hh) >= 1)

  block_sizes <- vapply(hh, length, integer(1))
  all_ids     <- sort(unique(unlist(hh)))
  n_total     <- length(all_ids)
  target      <- round(prop * n_total)

  best_sel   <- rep(FALSE, length(hh))
  best_sum   <- 0L
  best_diff  <- Inf

  for (it in seq_len(iterations)) {
    ord <- sample(seq_along(hh))
    sel <- rep(FALSE, length(hh))
    s   <- 0L

    # Greedy: traverse blocks in random order; add a block if it makes the total closer to the target
    for (j in ord) {
      new_s <- s + block_sizes[j]
      if (abs(new_s - target) <= abs(s - target)) {
        sel[j] <- TRUE
        s <- new_s
      }
      # Small optimization: exit early if already very close to the target
      if (abs(s - target) <= 1L) break
    }

    d <- abs(s - target)
    if (d < best_diff || (d == best_diff && s > best_sum)) {
      best_diff <- d
      best_sum  <- s
      best_sel  <- sel
      if (best_diff == 0L) {
        # Exactly hit the target, exit early
        break
      }
    }
  }

  train_ids <- sort(unique(unlist(hh[best_sel])))
  test_ids  <- setdiff(all_ids, train_ids)

  # Generate split vector sorted by sample IDs
  split <- rep(NA_character_, max(all_ids))
  split[train_ids] <- "train"
  split[test_ids]  <- "test"
  names(split) <- seq_along(split)

  list(
    train_ids = train_ids,
    test_ids  = test_ids,
    split     = split,
    sizes     = c(train = length(train_ids), test = length(test_ids), total = n_total),
    target    = c(target_train = target, achieved_train = length(train_ids), diff = length(train_ids) - target),
    block_sizes = block_sizes
  )
}

#' Merge a List of Results by Component
#'
#' This function takes a list of results (where each result is a named list or
#' vector with the same structure) and merges them by concatenating
#' corresponding components.
#'
#' @param res_list A list of named lists or vectors to be merged. All elements
#'   should have the same names/structure.
#'
#' @return A named list where each element contains the concatenated data
#'   from the corresponding components of the input list.
#' @export
merge_results <- function(res_list) {
  # Get names from the first element of the list
  nm <- names(res_list[[1]])

  # For each name, extract and concatenate values across all elements in res_list
  setNames(
    lapply(nm, function(x) do.call(c, lapply(res_list, `[[`, x))),
    nm
  )
}


#' Fix Constant Columns in a Matrix
#'
#' This function identifies columns that are constant within a specified subset
#' of rows and introduces artificial variation into those columns. This is often
#' used to prevent numerical issues (like singular matrices) in downstream
#' statistical tests.
#'
#' @param g_matrix A numeric matrix (e.g., a genotype matrix).
#' @param obs_idx An integer vector of row indices representing the observed samples to check.
#'
#' @return The modified \code{g_matrix} where constant columns in the observed
#'   subset have been perturbed.
#' @export
fix_constant_columns <- function(g_matrix, obs_idx) {
  # Extract the observed subset of the matrix
  g_obs <- as.matrix(g_matrix[obs_idx, ])

  # Identify columns where all values are identical to the first row
  ref_row <- g_obs[1, , drop = FALSE]
  is_const <- colSums(g_obs == matrix(rep(ref_row, each = nrow(g_obs)), nrow = nrow(g_obs))) == nrow(g_obs)

  # If constant columns exist, introduce variation
  if (any(is_const)) {
    const_cols <- which(is_const)

    # For each constant column, randomly pick 6 rows from the observed indices
    sampled_rows <- replicate(length(const_cols), sample(obs_idx, 6), simplify = FALSE)

    # Create an index matrix for assignment
    indices <- cbind(rows = unlist(sampled_rows), cols = rep(const_cols, each = 6))

    # Assign a pattern of values to break the constancy
    g_matrix[indices] <- rep(c(2, 1, 1, 1, 1, 2), length(const_cols))
  }

  return(g_matrix)
}


#' Rank-based Inverse Normal Transformation (INT)
#'
#' This function performs a rank-based inverse normal transformation (Blom's method by default)
#' on a specific phenotype column within a data frame. It handles missing values by
#' calculating the transformation based on the complete cases only.
#'
#' @param data A data frame containing the phenotype to be transformed.
#' @param pheno A character string specifying the name of the column in \code{data} to transform.
#' @param k A numeric constant used in the transformation formula (e.g., Blom's use 0.375). Default is 0.375.
#'
#' @return The original data frame with an additional column \code{int} containing the transformed values.
#' @export
INT <- function(data, pheno, k = 0.375){
  # Identify indices of missing values
  missing_index <- which(is.na(data[[pheno]]))

  if(length(missing_index) == 0){
    # All values are complete: use the whole data
    n <- nrow(data)
    r <- rank(data[[pheno]])
    data$int <- qnorm((r - k) / (n - 2 * k + 1))
  } else {
    # Some missing values exist: compute for complete cases only
    data.complete <- data[-missing_index, ]
    n <- nrow(data.complete)
    r <- rank(data.complete[[pheno]])
    data$int <- NA  # Initialize with NA
    data$int[-missing_index] <- qnorm((r - k) / (n - 2 * k + 1))
  }

  return(data)
}


#' Find Block Structures in a Sparse Matrix
#'
#' This function identifies independent block structures within a sparse matrix
#' (typically a precision or covariance matrix). It uses a vectorized approach
#' to determine the boundaries where the matrix can be partitioned into
#' non-overlapping sub-matrices.
#'
#' @param A A sparse matrix (usually of class \code{dgCMatrix} or similar).
#'
#' @return A list where each element is an integer vector containing the
#'   row/column indices for a specific block.
#' @export
find_blocks_vectorized <- function(A) {
  n <- nrow(A)
  # Extract triplet information (i, j, x) of non-zero elements
  sm <- summary(A)

  # Calculate the maximum non-zero column index for each row (0 if the row is all zeros)
  max_j <- rep(0, n)
  tmp <- tapply(sm$j, sm$i, max)
  max_j[as.integer(names(tmp))] <- tmp

  # For all-zero rows, set the value to its row index (indicating the row only affects itself)
  zero_rows <- which(max_j == 0)
  if(length(zero_rows) > 0){
    max_j[zero_rows] <- zero_rows
  }

  # Use cummax() to calculate the "furthest influence range" for each row
  cum_max <- cummax(max_j)

  # Identify block boundaries where the row index equals the cumulative maximum
  block_boundaries <- which(seq_len(n) == cum_max)

  # Partition into blocks based on identified boundaries
  blocks <- vector("list", length(block_boundaries))
  start <- 1
  for (i in seq_along(block_boundaries)) {
    end <- block_boundaries[i]
    blocks[[i]] <- start:end
    start <- end + 1
  }

  return(blocks)
}


#' Find Block Structures in a Sparse Matrix with a Threshold
#'
#' This function identifies independent block structures within a sparse matrix
#' (typically a precision or covariance matrix) after filtering out elements
#' whose absolute values are below a specified threshold. This helps in
#' partitioning the matrix into smaller, manageable sub-matrices.
#'
#' @param A A sparse matrix (usually of class \code{dgCMatrix} or similar).
#' @param threshold Numeric. Elements with absolute values smaller than this
#'   threshold are treated as zero. Default is 0.005.
#'
#' @return A list where each element is an integer vector containing the
#'   row/column indices for a specific block.
#' @export
find_blocks_vectorized_threshold <- function(A, threshold = 0.005) {
  n <- nrow(A)
  # Extract triplet information (i, j, x) and keep only elements >= threshold
  sm <- summary(A)
  sm <- sm[abs(sm$x) >= threshold, ]

  # Calculate the maximum non-zero column index for each row (0 if row is all zeros after filtering)
  max_j <- rep(0, n)
  if(nrow(sm) > 0) {
    tmp <- tapply(sm$j, sm$i, max)
    max_j[as.integer(names(tmp))] <- tmp
  }

  # For all-zero rows, set the value to its row index (indicating the row only affects itself)
  zero_rows <- which(max_j == 0)
  if (length(zero_rows) > 0) {
    max_j[zero_rows] <- zero_rows
  }

  # Calculate the "furthest influence range" for each row
  cum_max <- cummax(max_j)

  # Identify block boundaries where the row index equals the cumulative maximum
  block_boundaries <- which(seq_len(n) == cum_max)

  # Partition into blocks based on identified boundaries
  blocks <- vector("list", length(block_boundaries))
  start <- 1
  for (i in seq_along(block_boundaries)) {
    end <- block_boundaries[i]
    blocks[[i]] <- start:end
    start <- end + 1
  }

  return(blocks)
}


#' Construct a Block Matrix from Four Quadrants
#'
#' This function creates a large matrix by combining four individual matrices
#' representing the top-left, top-right, bottom-left, and bottom-right blocks.
#'
#' @param tl The top-left matrix block.
#' @param tr The top-right matrix block.
#' @param bl The bottom-left matrix block.
#' @param br The bottom-right matrix block.
#'
#' @return A single matrix formed by binding the four blocks together.
#' @export
block_matrix <- function(tl, tr, bl, br) {
  # Combine top blocks horizontally, combine bottom blocks horizontally,
  # and then bind the two resulting rows vertically.
  rbind(cbind(tl, tr), cbind(bl, br))
}


#' Parallel Block-wise Matrix Inversion
#'
#' This function performs an efficient inversion of a large sparse matrix by
#' partitioning it into independent blocks. It merges smaller blocks to
#' optimize parallel processing and uses Cholesky decomposition for
#' numerical stability.
#'
#' @param wait_matrix A large square matrix (typically a sparse \code{dgCMatrix}) to be inverted.
#' @param threshold Integer. The maximum number of rows/columns for a merged block
#'   to optimize memory and computation. Default is 1000.
#'
#' @return A sparse block-diagonal matrix representing the inverse of the input matrix.
#' @export
#' @import Matrix
#' @import parallel
matrix_inv_block <- function(wait_matrix, threshold = 1000){
  # Identify independent block structures within the matrix
  blocks <- find_blocks_vectorized(wait_matrix)

  merged_blocks <- list()
  current_block <- blocks[[1]]

  # Iterate through blocks and merge them if their combined size is below the threshold
  for (i in 2:length(blocks)) {
    if (length(current_block) + length(blocks[[i]]) < threshold) {
      current_block <- c(current_block, blocks[[i]])
    } else {
      # Otherwise, add the current block to the merged list and start a new one
      merged_blocks <- c(merged_blocks, list(current_block))
      current_block <- blocks[[i]]
    }
  }
  # Append the last accumulated block
  merged_blocks <- c(merged_blocks, list(current_block))

  # Set the number of cores for parallel computation
  ncores <- 2

  # Perform inversion for each block in parallel using Cholesky decomposition
  block_inv_list <- mclapply(merged_blocks, function(block) {
    subSigma <- wait_matrix[block, block]
    chol_subSigma <- chol(subSigma)
    chol2inv(chol_subSigma)
  }, mc.cores = ncores)

  # Combine individual block inverses into a single sparse block-diagonal matrix
  inv_wait_matrix <- bdiag(block_inv_list)

  return(inv_wait_matrix)
}

#' Invert a Relationship Matrix with Adaptive Block Detection
#'
#' This function performs block-wise inversion of a relationship matrix (A-matrix).
#' It features an adaptive thresholding mechanism that increases the sparseness
#' threshold until the maximum block size is manageable (<= 1000) for inversion.
#'
#' @param Amatrix A large square matrix (typically a sparse genetic relationship matrix).
#' @param thr Numeric. The initial threshold for detecting independent blocks. Default is 0.006.
#' @param threshold Integer. The maximum number of rows/columns for merged blocks
#'   during parallel computation. Default is 1000.
#'
#' @return A sparse block-diagonal matrix representing the inverse of the (possibly thresholded) input matrix.
#' @export
#' @import Matrix
#' @import parallel
matrix_inv_Amatrix <- function(Amatrix, thr = 0.006, threshold = 1000){
  # Detect independent blocks using an initial threshold
  blocks_A <- find_blocks_vectorized_threshold(Amatrix, threshold = thr)

  # Adaptive thresholding: increase thr if the largest block is too big for efficient inversion
  ja <- max(sapply(blocks_A, length))
  while(ja > 1000){
    thr <- thr + 0.001
    blocks_A <- find_blocks_vectorized_threshold(Amatrix, threshold = thr)
    ja <- max(sapply(blocks_A, length))
  }

  merged_blocks <- list()
  current_block <- blocks_A[[1]]

  # Merge adjacent blocks to reach the computational threshold size
  if(length(blocks_A) > 1 ){
    for (i in 2:length(blocks_A)) {
      if (length(current_block) + length(blocks_A[[i]]) < threshold) {
        current_block <- c(current_block, blocks_A[[i]])
      } else {
        # Otherwise, add the current block to the merged list and start a new one
        merged_blocks <- c(merged_blocks, list(current_block))
        current_block <- blocks_A[[i]]
      }
    }
  }

  # Append the last accumulated block
  merged_blocks <- c(merged_blocks, list(current_block))

  # Parallel computation using 2 cores
  ncores <- 2
  block_inv_list <- mclapply(merged_blocks, function(block) {
    subA <- Amatrix[block, block]
    # Solve the sub-matrix (Standard matrix inversion)
    solve(subA)
  }, mc.cores = ncores)

  # Combine block inverses into a single sparse block-diagonal matrix
  V11 <- bdiag(block_inv_list)

  return(V11)
}

#' Compute Score Test Statistic and P-value
#'
#' This function calculates the chi-square score statistic and its corresponding
#' negative log10 p-value based on the score vector and its variance.
#'
#' @param Ug Numeric. The score vector (numerator of the score test).
#' @param Vu Numeric. The variance of the score vector (denominator of the score test).
#'
#' @return A list containing:
#' \itemize{
#'   \item{score}: The computed chi-square score statistic.
#'   \item{negative_log10_pval}: The p-value on a -log10 scale.
#' }
#' @export
compute_score <- function(Ug, Vu) {
  # Calculate the chi-square statistic with 1 degree of freedom
  score <- as.numeric(Ug^2 / Vu)

  # Compute the -log10 p-value using the chi-square distribution
  # log.p = TRUE and division by log(10) is used for better numerical stability
  negative_log10_pval <- -pchisq(score, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)

  list(score = score, negative_log10_pval = negative_log10_pval)
}

#' Data Generating Process (DGP) for Syn-PALM Simulation
#'
#' This function generates simulated proteomic data (observed and synthetic) for
#' Evaluating Syn-PALM. It creates a population with cryptic relatedness using a
#' block-wise Genetic Relatedness Matrix (GRM), where individuals are grouped in
#' pairs with a kinship coefficient of 0.5 (simulating siblings or similar pairings).
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the proteomic data.
#' @param rho Correlation coefficient between the target protein (Y) and the surrogate/synthetic protein (S).
#' @param maf Minor Allele Frequency for the simulated SNP. Default is 0.25.
#' @param tauT Variance component of the random effect for the target protein.
#' @param tauS Variance component of the random effect for the surrogate protein.
#' @param sigmaT Variance component of the idiosyncratic error for the target protein.
#' @param sigmaS Variance component of the idiosyncratic error for the surrogate protein.
#' @param pve_g Proportion of Variance Explained by the genotype (heritability). Default is 0 (for Type I error control).
#' @param pve_x Proportion of Variance Explained by the covariate X. Default is 0.10.
#'
#' @return A list containing:
#' \itemize{
#'   \item{G_all}: Genotype vector for all samples.
#'   \item{X_all}: Covariate vector for all samples.
#'   \item{S}: Synthetic/Surrogate phenotype vector for all samples.
#'   \item{Y}: Oracle (fully observed) phenotype vector for all samples.
#'   \item{Y_obs}: Partially observed phenotype vector (NA for missing values).
#'   \item{GRM}: The sparse block-wise Genetic Relatedness Matrix.
#' }
#' @export
#' @import Matrix
DGP <- function(n_obs, miss, rho, maf = 0.25, tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,
                pve_g = 0, pve_x = 0.10) {

  ## Derive cross-trait variance components based on rho
  # Decompose the correlation into random effect (tauTS) and error (sigmaTS) components
  tauTS <- sqrt((tauT^2 + sigmaT^2) * (tauS^2 + sigmaS^2)) * rho * 3/5
  sigmaTS <- sqrt((tauT^2 + sigmaT^2) * (tauS^2 + sigmaS^2)) * rho * 2/5

  # Determine total sample size based on the target number of observed samples
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Generate Genotype (G) and Covariate (X)
  G_all <- stats::rbinom(n = n_all, size = 2, prob = maf)
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all, X_all)

  # Calculate total variance components for scaling fixed effects
  r <- tauT^2 + sigmaT^2
  V_total <- r / (1 - pve_g - pve_x)

  # Set fixed effect coefficients (Beta) based on PVE
  beta_G <- sqrt(pve_g * V_total / var(G_all))
  beta_X <- sqrt(pve_x * V_total / var(X_all))

  # Generate a block-wise GRM matrix (pairs of individuals with 0.5 relatedness)
  blocks <- seq(1, n_all - 1, by = 2)
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)
  values <- rep(c(1, 0.5, 0.5, 1), times = length(blocks))

  GRM <- sparseMatrix(i = rows, j = cols, x = values, dims = c(n_all, n_all))

  GRM_obs <- GRM[1:n_obs, 1:n_obs]
  GRM_obs_all <- GRM[1:n_obs, ]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all

  # Fixed effect component of the mean
  mu_all <- Z_all %*% c(beta_G, beta_X)

  ## Generate observed Y and surrogate S using the Joint Distribution
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)

  I_matrix <- sparseMatrix(
    i = 1:n_obs, j = 1:n_obs, x = sigmaTS,
    dims = c(n_obs, n_all)
  )

  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)

  # Cholesky decomposition for Sigma22 inversion
  chol_Sigma22 <- chol(Sigma22)
  inv_Sigma22 <- chol2inv(chol_Sigma22)

  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22

  # Generate Surrogate Phenotype S
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))

  # Generate Observed Phenotype Y_obs conditional on S
  mu_cond <- mu_all[1:n_obs] + Sigma12 %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond <- Sigma11 - Sigma12 %*% inv_Sigma22 %*% Sigma21

  L_cond <- chol(Sigma_cond)
  epsl <- rnorm(n_all)

  Y_obs <- rep(NA, n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])

  ### Generate Oracle Phenotype Y (Full observation benchmark)
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )

  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)

  mu_cond_oracle <- mu_all + Sigma12_oracle %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle %*% inv_Sigma22 %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)

  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)

  # Output results
  out <- list(
    G_all = G_all,
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )
  return(out)
}

#' Generate a Block-wise Genetic Relatedness Matrix (GRM)
#'
#' This helper function constructs a sparse block-diagonal GRM to simulate
#' cryptic relatedness in a study population. The matrix is composed of 2x2
#' blocks, where each block represents a pair of individuals with a specified
#' kinship correlation coefficient.
#'
#' @param n_obs Number of individuals with observed data.
#' @param miss Proportion of missingness, used to determine the total sample size.
#' @param GRMr Numeric value representing the relatedness coefficient between
#'   paired individuals (e.g., 0.5 for full siblings).
#'
#' @return A sparse square matrix of dimension \code{n_all x n_all} (where
#'   \code{n_all = n_obs / (1 - miss)}). The diagonal elements are 1, and
#'   off-diagonal elements within each 2x2 block are \code{GRMr}.
#' @export
#' @import Matrix
DGP_GRMr <- function(n_obs, miss, GRMr) {

  # Calculate total sample size based on observed count and missingness proportion
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Define the indices for block-wise structure (pairs of individuals)
  blocks <- seq(1, n_all - 1, by = 2)

  # Construct row and column indices for each 2x2 block
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)

  # Set the values: 1 for self-relatedness (diagonal), GRMr for paired-relatedness
  values <- rep(c(1, GRMr, GRMr, 1), times = length(blocks))

  # Create the sparse matrix efficiently using the Matrix package
  GRM <- sparseMatrix(i = rows, j = cols, x = values, dims = c(n_all, n_all))

  return(GRM)
}


#' Pre-calculate Variance-Covariance Components for Simulation Step 1
#'
#' This helper function pre-computes the required variance-covariance matrices
#' and their inversions for both the observed and synthetic proteomic data
#' components. It is designed to streamline the simulation process by
#' calculating static components once before iterating through multiple SNPs.
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the dataset.
#' @param tauT Variance component of the random effect for the target protein (Y).
#' @param tauS Variance component of the random effect for the surrogate protein (S).
#' @param sigmaT Variance component of the idiosyncratic error for the target protein.
#' @param sigmaS Variance component of the idiosyncratic error for the surrogate protein.
#' @param GRM A sparse Genetic Relatedness Matrix for all individuals in the study.
#'
#' @return A list containing pre-calculated matrix components:
#' \itemize{
#'   \item{Sigma11}: Covariance matrix for the observed samples.
#'   \item{Sigma22}: Covariance matrix for all samples (synthetic component).
#'   \item{inv_Sigma22}: Inverse of the \code{Sigma22} matrix.
#'   \item{Sigma11_oracle}: Full covariance matrix assuming all samples are observed.
#'   \item{chol_Sigma22}: Cholesky decomposition of the \code{Sigma22} matrix.
#' }
#' @export
#' @import Matrix
DGP_GRMr_nullmodel_step1 <- function(n_obs, miss, tauT = 0.7, tauS = 0.4,
                                     sigmaT = 0.7, sigmaS = 0.5, GRM) {

  # Determine total sample size
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Subset the GRM for observed samples
  GRM_obs <- GRM[1:n_obs, 1:n_obs]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs

  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all

  ## Pre-calculate covariance matrices to avoid redundant computation
  # Sigma11: Var(Y_obs)
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)

  # Sigma22: Var(S_all)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)

  # Efficient inversion of Sigma22 using Cholesky decomposition
  chol_Sigma22 <- chol(Sigma22)
  inv_Sigma22 <- chol2inv(chol_Sigma22)

  # Oracle covariance matrix (assuming full observation of Y)
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)

  # Collect results into a list for export
  params <- c("Sigma11", "Sigma22", "inv_Sigma22", "Sigma11_oracle", "chol_Sigma22")
  pre_vars_step1 <- mget(params, envir = environment())

  return(pre_vars_step1)
}


#' Pre-calculate Conditional Variance-Covariance Components for Simulation Step 2
#'
#' This function computes the cross-trait covariance matrices and the conditional
#' distribution parameters required for simulating phenotypes. It builds upon
#' the static components from Step 1, incorporating the correlation coefficient
#' (\code{rho}) to derive conditional means and variances for both the observed
#' and oracle scenarios.
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the dataset.
#' @param tauT Variance component of the random effect for the target protein.
#' @param tauS Variance component of the random effect for the surrogate protein.
#' @param sigmaT Variance component of the idiosyncratic error for the target protein.
#' @param sigmaS Variance component of the idiosyncratic error for the surrogate protein.
#' @param rho Correlation coefficient between the target and surrogate phenotypes.
#' @param GRM A sparse Genetic Relatedness Matrix for all samples.
#' @param pre_vars_step1 A list of pre-calculated variables from
#'   \code{DGP_GRMr_nullmodel_step1}.
#'
#' @return A list containing pre-calculated components for phenotype generation:
#' \itemize{
#'   \item{chol_Sigma22}: Cholesky decomposition of the surrogate covariance matrix.
#'   \item{Sigma12_Sigma22inv}: The product of cross-covariance and the inverse of
#'     surrogate covariance.
#'   \item{L_cond}: Cholesky decomposition of the conditional variance for
#'     observed samples.
#'   \item{Sigma12_oracle_Sigma22inv}: Oracle version of the projection matrix.
#'   \item{L_cond_oracle}: Cholesky decomposition of the oracle conditional variance.
#' }
#' @export
#' @import Matrix
DGP_GRMr_nullmodel_step2 <- function(n_obs, miss, tauT = 0.7, tauS = 0.4,
                                     sigmaT = 0.7, sigmaS = 0.5, rho, GRM,
                                     pre_vars_step1) {

  # Load static variables from Step 1 into the current environment
  list2env(pre_vars_step1, envir = environment())

  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Define subset indices for consistency
  GRM_obs_all <- GRM[1:n_obs, ]

  ## Generate cross-trait variance components based on rho
  # Allocation: 60% of correlation attributed to random effects, 40% to residuals
  tauTS <- sqrt((tauT^2 + sigmaT^2) * (tauS^2 + sigmaS^2)) * rho * 3/5
  sigmaTS <- sqrt((tauT^2 + sigmaT^2) * (tauS^2 + sigmaS^2)) * rho * 2/5

  # Construct the sparse cross-trait error covariance matrix
  I_matrix <- sparseMatrix(
    i = 1:n_obs, j = 1:n_obs, x = sigmaTS,
    dims = c(n_obs, n_all)
  )

  # Sigma12: Covariance between observed target and all surrogate phenotypes
  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)

  # Pre-calculate the projection matrix for the conditional mean
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22

  # Calculate the conditional variance matrix (Schur Complement) and its Cholesky factor
  Sigma_cond <- Sigma11 - Sigma12_Sigma22inv %*% Sigma21
  L_cond <- chol(Sigma_cond)

  ## Oracle calculations (assuming all target phenotypes are observed)
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )

  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)

  # Oracle projection matrix and conditional variance Cholesky factor
  Sigma12_oracle_Sigma22inv <- Sigma12_oracle %*% inv_Sigma22
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle_Sigma22inv %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)

  # Collect results for subsequent data generation (DGP)
  params <- c("chol_Sigma22", "Sigma12_Sigma22inv", "L_cond",
              "Sigma12_oracle_Sigma22inv", "L_cond_oracle")
  pre_vars_step2 <- mget(params, envir = environment())

  return(pre_vars_step2)
}

#' Generate Null Model Phenotypes for Simulation Step 3
#'
#' This function generates simulated phenotypes under the null hypothesis (no
#' genetic effect) for evaluating Type I error control. It utilizes the
#' pre-calculated variance-covariance decompositions from Step 2 to efficiently
#' produce surrogate and target phenotypes with the desired correlation and
#' relatedness structures.
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the dataset.
#' @param tauT Variance component of the random effect for the target protein.
#' @param sigmaT Variance component of the idiosyncratic error for the target protein.
#' @param pve_x Proportion of Variance Explained by the covariate X. Default is 0.10.
#' @param pre_vars_step2 A list of pre-calculated variables from
#'   \code{DGP_GRMr_nullmodel_step2} (e.g., \code{chol_Sigma22}, \code{L_cond}).
#' @param GRM A sparse Genetic Relatedness Matrix for all samples.
#'
#' @return A list containing the null-model simulated data:
#' \itemize{
#'   \item{X_all}: Covariate vector for all samples.
#'   \item{S}: Synthetic/Surrogate phenotype vector for all samples.
#'   \item{Y}: Oracle (fully observed) target phenotype vector for all samples.
#'   \item{Y_obs}: Partially observed target phenotype vector (NA for missing).
#'   \item{GRM}: The Genetic Relatedness Matrix used for simulation.
#' }
#' @export
DGP_GRMr_nullmodel_step3 <- function(n_obs, miss, tauT = 0.7, sigmaT = 0.7, pve_x = 0.10,
                                     pre_vars_step2, GRM) {

  # Load pre-calculated decomposition components into the environment
  list2env(pre_vars_step2, envir = environment())

  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Generate covariate X (standard normal) and set genotype G to zero (Null Model)
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(rep(0, n_all), X_all)

  # Calculate total phenotype variation to scale the fixed effect of X
  r <- tauT^2 + sigmaT^2
  V_total <- r / (1 - 0 - pve_x)

  # Set fixed effect coefficients: beta_G is 0 under the null
  beta_G <- 0
  beta_X <- sqrt(pve_x * V_total / var(X_all))

  # Compute the fixed effect component of the mean
  mu_all <- Z_all %*% c(beta_G, beta_X)

  ## Generate Surrogate Phenotype S using pre-calculated Cholesky factor
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))

  ## Generate Observed Target Phenotype Y_obs conditional on S
  # mu_cond follows the property of conditional multivariate normal distribution
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)

  Y_obs <- rep(NA, n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])

  ### Generate Oracle Target Phenotype Y (benchmark with no missingness)
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)

  # Consolidate results for GWAS evaluation
  out <- list(
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )
  return(out)
}

#' Generate Genotype Matrix for Type I Error Simulation
#'
#' This helper function generates a matrix of simulated genotypes for a specified
#' number of SNPs (chunk size) across all individuals. The genotypes are
#' simulated under the null hypothesis (independent of the phenotype) using a
#' binomial distribution based on a given Minor Allele Frequency (MAF).
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the dataset, used to calculate
#'   the total sample size.
#' @param chunk_size Integer; the number of SNPs to generate in this matrix
#'   (e.g., 100 or 1000).
#' @param maf Minor Allele Frequency (MAF) used for the binomial distribution
#'   (\code{size = 2}).
#'
#' @return A numeric matrix of dimensions \code{n_all x chunk_size}, where
#'   \code{n_all} is the total number of individuals (observed + missing).
#'   Each element represents the number of risk alleles (0, 1, or 2).
#' @export
typeIerror_generate_gmatrix <- function(n_obs, miss, chunk_size, maf) {

  # Calculate the number of missing samples and the resulting total sample size
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Generate the genotype matrix using a binomial distribution
  # size = 2 represents a diploid organism
  g_matrix <- matrix(rbinom(n_all * chunk_size, size = 2, prob = maf),
                     nrow = n_all,
                     ncol = chunk_size)

  return(g_matrix)
}

#' Generate Genotype Vector with Specific Kinship Structures
#'
#' This function simulates the inheritance process to generate genotypes for
#' pairs of individuals under various kinship settings. It uses a gamete-drawing
#' mechanism to ensure the resulting genetic relatedness (e.g., kinship
#' coefficient) matches biological expectations.
#'
#' @param n_obs Number of individuals with observed proteomic data.
#' @param miss Proportion of missingness in the dataset.
#' @param maf Minor Allele Frequency (MAF) for the simulated SNP.
#' @param kinsetting Character string specifying the kinship structure.
#'   Options include:
#'   \itemize{
#'     \item{\code{"unkin"}}: Unrelated individuals (standard binomial draw).
#'     \item{\code{"siblings"}}: Full siblings (sharing the same two parents).
#'     \item{\code{"cousin"}}: First cousins (parents are siblings).
#'     \item{\code{"halfsiblings"}}: Half-siblings (sharing one parent).
#'   }
#'
#' @return A numeric vector of genotypes (0, 1, or 2) of length \code{n_all},
#'   where individuals are interleaved to form pairs according to the
#'   specified \code{kinsetting}.
#' @export
generate_gmatrix_kinship <- function(n_obs, miss, maf, kinsetting) {

  # Calculate total sample size
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Helper function to draw a single gamete (0 or 1) from a diploid genotype (0, 1, or 2)
  draw_gamete <- function(g) rbinom(length(g), size = 1, prob = g / 2)

  # Generate founder genotypes (Grandparents/Parents level)
  g_grand_father <- rbinom(n_all/2, size = 2, prob = maf)
  g_grand_mother <- rbinom(n_all/2, size = 2, prob = maf)

  # Case 1: Unrelated individuals
  if (kinsetting == "unkin") {
    geno_vec <- c(rbind(g_grand_father, g_grand_mother))
  }

  # Case 2: Full siblings (Shared father and mother)
  if (kinsetting == "siblings") {
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    geno_vec <- c(rbind(g_father, g_uncle))
  }

  # Case 3: First cousins (Parents are siblings)
  if (kinsetting == "cousin") {
    # Two parents are siblings
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)

    # Marry unrelated founders
    g_mother <- rbinom(n_all/2, size = 2, prob = maf)
    g_me <- draw_gamete(g_father) + draw_gamete(g_mother)

    g_rand <- rbinom(n_all/2, size = 2, prob = maf)
    g_far_sister <- draw_gamete(g_uncle) + draw_gamete(g_rand)

    geno_vec <- c(rbind(g_me, g_far_sister))
  }

  # Case 4: Half-siblings (Shared father, different mothers)
  if (kinsetting == "halfsiblings") {
    # Child 1 from Parent A + Parent B
    g_uncle1 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)

    # Child 2 from Parent A + Parent C (new unrelated founder)
    g_grand_mother2 <- rbinom(n_all/2, size = 2, prob = maf)
    g_uncle2 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother2)

    geno_vec <- c(rbind(g_uncle1, g_uncle2))
  }

  return(geno_vec)
}

#' Data Generation Process (DGP) for Power Comparison (Mixed Model)
#'
#' This function generates phenotypes (Target Y and Surrogate S) conditioned on
#' a Genetic Relatedness Matrix (GRM) and specified effect sizes (PVE). It uses
#' pre-calculated Cholesky decompositions from Step 2 to efficiently generate
#' correlated traits, ensuring the simulated power reflects the underlying
#' genetic architecture.
#'
#' @param n_obs Number of individuals with observed phenotypes.
#' @param miss Proportion of missingness for the target phenotype.
#' @param tauT Variance component for the random effect in the target trait.
#' @param sigmaT Residual variance component for the target trait.
#' @param pve_x Proportion of variance explained by covariates.
#' @param pve_g Proportion of variance explained by the genetic variant (Effect Size).
#' @param maf Minor Allele Frequency for the simulated SNP.
#' @param kinsetting Kinship setting for genotype generation.
#' @param Genotype_each Optional: provide a pre-existing genotype vector.
#' @param pre_vars_step2 A list containing pre-calculated covariance matrices
#'   (e.g., chol_Sigma22, L_cond, Sigma12_Sigma22inv) from Step 2.
#' @param GRM The full Genetic Relatedness Matrix.
#'
#' @return A list containing the simulated datasets:
#' \itemize{
#'   \item{\code{G_all}}: Genotype vector for all samples.
#'   \item{\code{X_all}}: Covariate vector for all samples.
#'   \item{\code{S}}: Synthetic/Surrogate phenotypes for all samples.
#'   \item{\code{Y}}: Oracle (fully observed) target phenotypes.
#'   \item{\code{Y_obs}}: Partially observed target phenotypes (with NAs).
#'   \item{\code{GRM}}: The Genetic Relatedness Matrix used.
#' }
#' @export
DGP_GRMr_gmodel_step3 <- function(n_obs, miss, tauT = 0.7, sigmaT = 0.7, pve_x = 0, pve_g = 0.0001,
                                  maf, kinsetting, Genotype_each = NULL, pre_vars_step2, GRM) {

  # Unpack pre-calculated covariance structures into the current environment
  list2env(pre_vars_step2, envir = environment())

  # Calculate total sample size based on missingness rate
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # 1. Genotype and Covariate Generation
  if(!is.na(maf)){
    G_all <- generate_gmatrix_kinship(n_obs, miss, maf, kinsetting)
  }
  if(!is.null(Genotype_each)){
    G_all <- Genotype_each
  }

  X_all <- stats::rnorm(n = n_all)
  # Design matrix for fixed effects (SNP + Covariate)
  Z_all <- cbind(G_all, X_all)

  # 2. Variance Scaling
  # Calculate total phenotype variation based on desired PVE (Proportion of Variance Explained)
  r <- tauT^2 + sigmaT^2
  V_total <- r / (1 - pve_g - pve_x)

  # Derive fixed effect coefficients to match the specified PVE
  beta_G <- sqrt(pve_g * V_total / var(G_all))
  beta_X <- sqrt(pve_x * V_total / var(X_all))

  # Calculate the linear predictor (Fixed effects part)
  mu_all <- Z_all %*% c(beta_G, beta_X)

  # 3. Trait Generation via Conditional Gaussian Logic
  # Step A: Generate Surrogate Phenotype S
  # S ~ N(mu_all, Sigma22)
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))

  # Step B: Generate Observed Target Y conditioned on S
  # Y_obs | S ~ N(mu_cond, Sigma_cond)
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA, n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])

  # Step C: Generate Oracle Y (Full cohort) conditioned on S
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)

  # Package output
  out <- list(
    G_all = G_all,
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )

  return(out)
}

#' Generate Correlated Genotypes for Nuclear Families (Latent Variable Method)
#'
#' This function simulates a genotype matrix for a cohort composed of nuclear
#' families (pairs of related individuals). It uses a latent normal distribution
#' to introduce correlation (rho) between family members and then transforms
#' these values into discrete genotypes (0, 1, 2) based on the specified
#' Minor Allele Frequency (MAF).
#'
#' @param n_obs Number of individuals with observed phenotypes (used to derive total N).
#' @param miss Proportion of missingness to determine the total cohort size.
#' @param n_snps Number of SNPs to generate.
#' @param maf Minor Allele Frequency (default is 0.25).
#' @param rho Genetic correlation between family members (e.g., 0.5 for siblings).
#'
#' @return A standardized genotype matrix of size \eqn{N \times n\_snps}, where
#'   rows represent individuals and columns represent SNPs. The genotypes are
#'   centered and scaled to have mean 0 and variance 1.
#'
#' @details
#' The simulation follows these steps:
#' \enumerate{
#'   \item Generate independent standard normal variables for family members.
#'   \item Introduce correlation using the formula: \eqn{Z_{2} = \rho Z_{1} + \sqrt{1-\rho^2} Z_{new}}.
#'   \item Apply symmetric thresholds based on \code{qnorm(1 - maf)} to convert
#'         latent variables into 0, 1, and 2.
#'   \item Standardize the resulting matrix using the theoretical variance \eqn{2 \cdot MAF \cdot (1 - MAF)}.
#' }
#' @export
generate_genotypes_rho <- function(n_obs, miss, n_snps, maf = 0.25, rho = 0.5) {

  # Calculate total cohort size
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # Ensure an even number for pair-based family simulation
  n_families <- n_all / 2

  # 1. Generate independent latent normal variables
  # Z is a matrix where each row is a SNP and columns are independent samples
  Z <- matrix(rnorm(n_snps * 2 * n_families), nrow = n_snps)

  # 2. Induce correlation between family members
  # For each family (pair), Person 1 is independent; Person 2 is correlated to Person 1.
  Z_family <- matrix(0, nrow = n_snps, ncol = 2 * n_families)

  for (i in 1:n_families) {
    # First family member (independent)
    Z_family[, 2*i - 1] <- Z[, 2*i - 1]

    # Second family member (correlated with rho)
    Z_family[, 2*i] <- rho * Z[, 2*i - 1] + sqrt(1 - rho^2) * Z[, 2*i]
  }

  # 3. Transform latent variables to discrete genotypes (0, 1, 2)
  # Determine thresholds based on MAF under normality assumption
  threshold <- qnorm(1 - maf)
  genotypes <- matrix(0, nrow = 2 * n_families, ncol = n_snps)

  for (i in 1:n_snps) {
    # Apply symmetric bilateral thresholds:
    # 0: < -threshold (Major allele homozygote)
    # 2: > threshold  (Minor allele homozygote)
    # 1: between     (Heterozygote)
    upper_threshold <- threshold
    lower_threshold <- -threshold

    gt <- ifelse(Z_family[i,] < lower_threshold, 0,
                 ifelse(Z_family[i,] > upper_threshold, 2, 1))

    genotypes[, i] <- gt
  }

  # 4. Standardize genotypes (Mean = 0, Variance = 1)
  # Centering
  genotypes_centered <- scale(genotypes, center = TRUE, scale = FALSE)

  # Scaling using theoretical binomial variance
  genotypes_scaled <- genotypes_centered / sqrt(2 * maf * (1 - maf))

  # Final output: individuals as rows, SNPs as columns
  return(genotypes_scaled)
}

#' Data Generation Process (DGP) for Correlated Samples (Mixed Model)
#'
#' This function generates phenotypes (Target Y and Surrogate S) conditioned on
#' a pre-existing genotype matrix (e.g., simulated with specific family
#' correlations) and a Genetic Relatedness Matrix (GRM). It utilizes
#' pre-calculated conditional distributions to ensure that the simulated traits
#' maintain the intended variance components and genetic architecture.
#'
#' @param n_obs Number of individuals with observed phenotypes.
#' @param miss Proportion of missingness for the target phenotype.
#' @param tauT Variance component for the random effect in the target trait.
#' @param sigmaT Residual variance component for the target trait.
#' @param pve_x Proportion of variance explained by covariates.
#' @param pve_g Proportion of variance explained by the specific SNP (Effect Size).
#' @param Genotype_each A vector of genotypes for all samples (typically generated
#'   via \code{generate_genotypes_rho}).
#' @param pre_vars_step2 A list containing pre-calculated covariance matrices
#'   (e.g., \code{chol_Sigma22}, \code{L_cond}, \code{Sigma12_Sigma22inv}) from Step 2.
#' @param GRM The full Genetic Relatedness Matrix corresponding to the samples.
#'
#' @return A list containing the simulated dataset:
#' \itemize{
#'   \item{\code{G_all}}: The input genotype vector.
#'   \item{\code{X_all}}: Simulated independent normal covariates.
#'   \item{\code{S}}: Synthetic/Surrogate phenotypes for all samples.
#'   \item{\code{Y}}: Oracle (fully observed) target phenotypes.
#'   \item{\code{Y_obs}}: Partially observed target phenotypes (with NAs).
#'   \item{\code{GRM}}: The Genetic Relatedness Matrix used.
#' }
#' @export
DGP_GRMr_gmodel_step3_rho <- function(n_obs, miss, tauT = 0.7, sigmaT = 0.7, pve_x = 0, pve_g = 0.0001,
                                      Genotype_each, pre_vars_step2, GRM) {

  # Unpack pre-calculated covariance structures into the current environment
  list2env(pre_vars_step2, envir = environment())

  # Calculate total cohort size
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss

  # 1. Fixed Effects Construction
  # Use external genotypes (e.g., from family-based simulation)
  G_all <- Genotype_each
  X_all <- stats::rnorm(n = n_all)

  # Design matrix for fixed effects
  Z_all <- cbind(G_all, X_all)

  # 2. Variance and Coefficient Scaling
  # r represents the total variance contributed by the random effect and error
  r <- tauT^2 + sigmaT^2
  # Rescale to ensure fixed effects achieve the target PVE
  V_total <- r / (1 - pve_g - pve_x)

  # Derive coefficients based on the variance of the generated Genotypes/Covariates
  beta_G <- sqrt(pve_g * V_total / var(G_all))
  beta_X <- sqrt(pve_x * V_total / var(X_all))

  # Linear predictor (Fixed effect component)
  mu_all <- Z_all %*% c(beta_G, beta_X)

  # 3. Trait Generation via Conditional Sampling
  # Step A: Generate Surrogate Phenotype S based on pre-calculated Sigma22
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))

  # Step B: Generate Observed Target Y conditioned on S
  # Uses the Schur complement logic (Sigma12 * Sigma22^-1)
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA, n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])

  # Step C: Generate Oracle Y (Full cohort) conditioned on S
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)

  # Package output
  out <- list(
    G_all = G_all,
    X_all = X_all,
    S = S,
    Y = Y,
    Y_obs = Y_obs,
    GRM = GRM
  )

  return(out)
}

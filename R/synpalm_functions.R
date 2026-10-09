## ---------------------------------------------------------------------------
## SynPALM: function definitions.
##
## The code below is taken verbatim from the analysis source used to produce
## the manuscript results, so that the package and the published analysis are
## guaranteed to agree.
## ---------------------------------------------------------------------------

#' Clear and recreate a scratch directory
#'
#' Deletes \code{dir} and everything inside it if it exists, then creates it
#' again as an empty directory. Used to give each job a clean temporary
#' workspace before writing intermediate files.
#'
#' @param dir Character path to the directory to reset.
#'
#' @return Called for its side effect on the file system; returns \code{NULL}
#'   invisibly.
#' @export
reset_temp_dir <- function(dir) {
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE) # 删除现有目录
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE) # 重新创建
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
  # eid: 向量 (UKB ID 或其他 sample ID)
  # covariates: data.frame，行数 = length(eid)，列为附加协变量
  # file_out: 输出文件名
  
  n <- length(eid)
  
  # 基础列
  df <- data.frame(
    ID_1 = eid,
    ID_2 = eid,
    missing = 0
  )
  
  # 添加协变量
  if (!is.null(covariates)) {
    stopifnot(nrow(covariates) == n)
    df <- cbind(df, covariates)
  }
  
  # === 写文件 ===
  # 第一行：列名
  cat(paste(colnames(df), collapse = " "), "\n", file = file_out)
  
  # 第二行：类型说明
  #   ID_1, ID_2, missing -> "0"
  #   其他列 -> "D" (表示定性变量 / covariate)
  type_line <- c("0", "0", "0", rep("D", ncol(df) - 3))
  cat(paste(type_line, collapse = " "), "\n", file = file_out, append = TRUE)
  
  # 后续行：样本信息
  write.table(
    df,
    file = file_out,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
  
  message("Done: .sample file written to: ", file_out)
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
    
    # 贪心：按随机顺序遍历 block，若加入能让总数更接近 target 就加入
    for (j in ord) {
      new_s <- s + block_sizes[j]
      if (abs(new_s - target) <= abs(s - target)) {
        sel[j] <- TRUE
        s <- new_s
      }
      # 小优化：已经非常接近就可以提前结束
      if (abs(s - target) <= 1L) break
    }
    
    d <- abs(s - target)
    if (d < best_diff || (d == best_diff && s > best_sum)) {
      best_diff <- d
      best_sum  <- s
      best_sel  <- sel
      if (best_diff == 0L) {
        # 已经正好命中目标，直接退出
        break
      }
    }
  }
  
  train_ids <- sort(unique(unlist(hh[best_sel])))
  test_ids  <- setdiff(all_ids, train_ids)
  
  # 生成按样本编号排序的 split 向量
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
  nm <- names(res_list[[1]])
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
  g_obs <- as.matrix(g_matrix[obs_idx, ])
  ref_row <- g_obs[1, , drop = FALSE]
  is_const <- colSums(g_obs == matrix(rep(ref_row, each = nrow(g_obs)), nrow = nrow(g_obs))) == nrow(g_obs)
  
  if (any(is_const)) {
    const_cols <- which(is_const)
    sampled_rows <- replicate(length(const_cols), sample(obs_idx, 6), simplify = FALSE)
    indices <- cbind(rows = unlist(sampled_rows), cols = rep(const_cols, each = 6))
    g_matrix[indices] <- rep(c(2, 1, 1, 1, 1, 2), length(const_cols))
  }
  g_matrix
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
    data$int <- NA  # initialize
    data$int[-missing_index] <- qnorm((r - k) / (n - 2 * k + 1))
  }
  
  return(data)
}

#' Read a SynPALM numerical-stability option
#'
#' Looks up \code{getOption("synsurrg.<name>")}, falling back to \code{default}
#' when the option is unset or \code{NULL}.
#'
#' @param name Option name without the \code{"synsurrg."} prefix.
#' @param default Value returned when the option is not set.
#'
#' @return The option value, or \code{default}.
#' @noRd
.ss_opt <- function(name, default) {
  v <- getOption(paste0("synsurrg.", name), default)
  if (is.null(v)) default else v
}

#' Number of cores used for block-wise inversion
#'
#' Uses \code{getOption("synsurrg.ncores")} if set, otherwise the
#' \code{SLURM_CPUS_PER_TASK} environment variable, otherwise 2.
#'
#' @return A positive integer.
#' @noRd
.ss_cores <- function() {
  v <- getOption("synsurrg.ncores", NULL)
  if (!is.null(v) && is.finite(v) && v >= 1) return(as.integer(v))
  x <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (!is.na(x) && x >= 1L) return(x)
  2L
}

#' Whether to print variance-component and block-inversion diagnostics
#'
#' @return \code{TRUE} unless \code{options(synsurrg.verbose = FALSE)}.
#' @noRd
.ss_verbose <- function() isTRUE(.ss_opt("verbose", TRUE))

#' Whether non-PSD estimates should raise an error instead of being projected
#'
#' @return \code{TRUE} if \code{options(synsurrg.on_nonpsd = "stop")}.
#' @noRd
.ss_strict  <- function() identical(.ss_opt("on_nonpsd", "project"), "stop")

#' Project a Symmetric 2x2 Matrix onto the PSD Cone
#'
#' Truncates negative eigenvalues of \code{[a, ab; ab, b]} to zero, applies a
#' lower bound to the diagonal, and pulls the implied correlation just inside
#' the unit circle so the result is strictly positive definite. All three steps
#' are no-ops for a matrix that is already comfortably positive definite, so
#' well-behaved estimates are returned unchanged. Non-finite entries are set to
#' zero before projecting.
#'
#' @param a,b Diagonal entries.
#' @param ab Off-diagonal entry.
#' @param floor_diag Lower bound applied to the diagonal after the eigenvalue
#'   projection (adding a non-negative diagonal preserves PSD).
#' @param shrink Relative amount by which \code{|ab|} is kept below
#'   \code{sqrt(a * b)}.
#'
#' @return A list with the projected entries \code{a}, \code{ab}, \code{b};
#'   a logical \code{changed}; and the flags \code{truncated}, \code{floored}
#'   and \code{shrunk} recording which step modified the matrix.
#' @keywords internal
#' @export
psd_project_2x2 <- function(a, ab, b, floor_diag = 0, shrink = 1e-6) {
  M   <- matrix(c(a, ab, ab, b), 2L, 2L)
  bad <- !is.finite(M)
  M[bad] <- 0
  
  e <- eigen(M, symmetric = TRUE)
  truncated <- any(e$values < 0)
  if (truncated) {
    e$values[e$values < 0] <- 0
    M <- e$vectors %*% diag(e$values, 2L) %*% t(e$vectors)
  }
  
  floored <- any(diag(M) < floor_diag)
  if (floored) diag(M) <- pmax(diag(M), floor_diag)
  
  lim      <- (1 - shrink) * sqrt(max(M[1, 1] * M[2, 2], 0))
  shrunk   <- abs(M[1, 2]) > lim
  if (shrunk) M[1, 2] <- M[2, 1] <- sign(M[1, 2]) * lim
  
  list(a = M[1, 1], ab = M[1, 2], b = M[2, 2],
       changed = truncated || floored || shrunk || any(bad),
       truncated = truncated, floored = floored, shrunk = shrunk)
}

#' Constrain Bivariate Variance Components to a Valid Covariance
#'
#' The Haseman-Elston moment estimators of the target (T) and synthetic (S)
#' variance components are unconstrained. The joint covariance
#' \eqn{G \otimes GRM + E \otimes I} is positive definite iff the genetic
#' 2x2 matrix \eqn{G} is PSD and the residual 2x2 matrix \eqn{E} is PD, and
#' flooring the six numbers one at a time does not guarantee that. This
#' function projects \eqn{G} and \eqn{E} onto the PSD cone with
#' \code{\link{psd_project_2x2}}, flooring the residual variances at
#' \code{getOption("synsurrg.ridge_e", 0.01)} times the trait scale. This
#' guarantees that \code{Sigma11}, \code{Sigma22} and the Schur complement
#' \code{A} are all positive definite.
#'
#' With \code{options(synsurrg.on_nonpsd = "stop")} an error is raised instead
#' of projecting.
#'
#' @param tau_T2,tau_TS,tau_S2 Raw genetic variance and covariance estimates.
#' @param sigma_T2,sigma_TS,sigma_S2 Raw residual variance and covariance
#'   estimates.
#' @param var_T,var_S Mean squared residuals of the two traits, used to scale
#'   the residual-variance floor.
#' @param label Character label used in messages.
#'
#' @return A list with the six constrained components and \code{vc}, a list
#'   holding the \code{raw} and \code{used} values, the projection flags, and
#'   the implied genetic and residual correlations \code{rho_g}, \code{rho_e}.
#' @keywords internal
#' @export
constrain_vc_bivariate <- function(tau_T2, tau_TS, tau_S2,
                                   sigma_T2, sigma_TS, sigma_S2,
                                   var_T = 1, var_S = 1, label = "") {
  raw <- c(tau_T2 = tau_T2, tau_TS = tau_TS, tau_S2 = tau_S2,
           sigma_T2 = sigma_T2, sigma_TS = sigma_TS, sigma_S2 = sigma_S2)
  
  eps   <- .ss_opt("ridge_e", 0.01)
  scale <- mean(c(abs(var_T), abs(var_S)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  
  g <- psd_project_2x2(tau_T2,   tau_TS,   tau_S2,   floor_diag = 0)
  r <- psd_project_2x2(sigma_T2, sigma_TS, sigma_S2, floor_diag = eps * scale)
  
  ## "none": keep the raw estimates as the pre-2026-08-22 code did, flooring
  ## only sigma_T2 at 0.01. The projections above are still computed so that
  ## vc records whether they would have changed anything.
  if (identical(.ss_opt("on_nonpsd", "project"), "none")) {
    used <- raw
    used[["sigma_T2"]] <- max(used[["sigma_T2"]], 0.01)
    vc <- list(raw = raw, used = used, mode = "none",
               projected_G = FALSE, projected_E = FALSE,
               would_project_G = g$changed, would_project_E = r$changed,
               detail_G = g[c("truncated", "floored", "shrunk")],
               detail_E = r[c("truncated", "floored", "shrunk")],
               rho_g = used[["tau_TS"]]   / sqrt(max(used[["tau_T2"]]   * used[["tau_S2"]],   .Machine$double.eps)),
               rho_e = used[["sigma_TS"]] / sqrt(max(used[["sigma_T2"]] * used[["sigma_S2"]], .Machine$double.eps)),
               label = label)
    if (.ss_verbose() && (g$changed || r$changed))
      message(sprintf("[%s] VC not PSD (G %s, E %s); kept unconstrained (on_nonpsd = \"none\")",
                      label, if (g$changed) "invalid" else "ok",
                      if (r$changed) "invalid" else "ok"))
    return(list(tau_T2 = used[["tau_T2"]], tau_TS = used[["tau_TS"]], tau_S2 = used[["tau_S2"]],
                sigma_T2 = used[["sigma_T2"]], sigma_TS = used[["sigma_TS"]],
                sigma_S2 = used[["sigma_S2"]], vc = vc))
  }
  
  used <- c(tau_T2 = g$a, tau_TS = g$ab, tau_S2 = g$b,
            sigma_T2 = r$a, sigma_TS = r$ab, sigma_S2 = r$b)
  
  if ((g$changed || r$changed) && .ss_strict())
    stop(sprintf("[%s] variance components not PSD (G %s, E %s): %s",
                 label, if (g$changed) "adjusted" else "ok",
                 if (r$changed) "adjusted" else "ok",
                 paste(sprintf("%s=%+.4g", names(raw), raw), collapse = " ")),
         call. = FALSE)
  
  vc <- list(raw = raw, used = used,
             projected_G = g$changed, projected_E = r$changed,
             detail_G = g[c("truncated", "floored", "shrunk")],
             detail_E = r[c("truncated", "floored", "shrunk")],
             rho_g = used[["tau_TS"]]   / sqrt(max(used[["tau_T2"]]   * used[["tau_S2"]],   .Machine$double.eps)),
             rho_e = used[["sigma_TS"]] / sqrt(max(used[["sigma_T2"]] * used[["sigma_S2"]], .Machine$double.eps)),
             label = label)
  
  if (.ss_verbose() && (g$changed || r$changed))
    message(sprintf("[%s] VC projected onto PSD cone\n  raw : %s\n  used: %s\n  rho_g = %.3f, rho_e = %.3f",
                    label,
                    paste(sprintf("%s=%+.4g", names(raw),  raw),  collapse = " "),
                    paste(sprintf("%s=%+.4g", names(used), used), collapse = " "),
                    vc$rho_g, vc$rho_e))
  
  list(tau_T2 = used[["tau_T2"]], tau_TS = used[["tau_TS"]], tau_S2 = used[["tau_S2"]],
       sigma_T2 = used[["sigma_T2"]], sigma_TS = used[["sigma_TS"]],
       sigma_S2 = used[["sigma_S2"]], vc = vc)
}

#' Constrain Univariate Variance Components
#'
#' \eqn{\tau^2 GRM + \sigma^2 I} is positive definite iff \eqn{\tau^2 \ge 0}
#' and \eqn{\sigma^2 > 0}, so a simple floor is the correct constraint here.
#' \eqn{\sigma^2} is floored at \code{getOption("synsurrg.ridge_e", 0.01)}
#' times \code{var_T}.
#'
#' @param tau2 Raw genetic variance estimate.
#' @param sigma2 Raw residual variance estimate.
#' @param var_T Mean squared residual of the trait, used to scale the floor.
#' @param label Character label used in messages.
#'
#' @return A list with \code{tau2}, \code{sigma2} and \code{vc}, a list
#'   holding the \code{raw} and \code{used} values and the flooring flags.
#' @keywords internal
#' @export
constrain_vc_univariate <- function(tau2, sigma2, var_T = 1, label = "") {
  raw <- c(tau_T2 = tau2, sigma_T2 = sigma2)
  
  eps   <- .ss_opt("ridge_e", 0.01)
  scale <- abs(var_T)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  
  ## "none": pre-2026-08-22 behaviour, sigma2 floored at 0.01 and nothing else.
  if (identical(.ss_opt("on_nonpsd", "project"), "none")) {
    sigma2 <- max(sigma2, 0.01)
    return(list(tau2 = tau2, sigma2 = sigma2,
                vc = list(raw = raw, used = c(tau_T2 = tau2, sigma_T2 = sigma2), mode = "none",
                          floored_tau = FALSE, floored_sigma = sigma2 != raw[["sigma_T2"]],
                          label = label)))
  }
  
  bad_tau   <- !is.finite(tau2)   || tau2   < 0
  bad_sigma <- !is.finite(sigma2) || sigma2 < eps * scale
  tau2   <- if (bad_tau)   0             else tau2
  sigma2 <- if (bad_sigma) eps * scale   else sigma2
  
  if ((bad_tau || bad_sigma) && .ss_strict())
    stop(sprintf("[%s] variance components out of range: %s", label,
                 paste(sprintf("%s=%+.4g", names(raw), raw), collapse = " ")),
         call. = FALSE)
  
  if (.ss_verbose() && (bad_tau || bad_sigma))
    message(sprintf("[%s] VC floored: tau_T2 %+.4g -> %.4g, sigma_T2 %+.4g -> %.4g",
                    label, raw[["tau_T2"]], tau2, raw[["sigma_T2"]], sigma2))
  
  list(tau2 = tau2, sigma2 = sigma2,
       vc = list(raw = raw, used = c(tau_T2 = tau2, sigma_T2 = sigma2),
                 floored_tau = bad_tau, floored_sigma = bad_sigma, label = label))
}

#' Greedily merge consecutive blocks up to a size limit
#'
#' @param blocks List of index vectors from \code{find_blocks_vectorized()}.
#' @param threshold Merged blocks are kept strictly smaller than this size.
#'
#' @return A list of merged index vectors.
#' @noRd
.ss_merge_blocks <- function(blocks, threshold) {
  if (!length(blocks)) stop(".ss_merge_blocks: empty block list")
  merged <- list()
  cur <- blocks[[1L]]
  if (length(blocks) >= 2L) {          # 2:length(blocks) was c(2,1) when length==1
    for (i in 2:length(blocks)) {
      if (length(cur) + length(blocks[[i]]) < threshold) {
        cur <- c(cur, blocks[[i]])
      } else {
        merged <- c(merged, list(cur))
        cur <- blocks[[i]]
      }
    }
  }
  c(merged, list(cur))
}

#' Invert one block with escalating fallbacks
#'
#' Tries a Cholesky inverse, then Cholesky with an escalating ridge, then an
#' eigenvalue-floored inverse.
#'
#' @param M Square block to invert.
#' @param max_tries Number of Cholesky attempts (the first without a ridge).
#' @param dense_max Blocks up to this size are inverted as dense matrices.
#' @param diag_eigen_max Blocks up to this size get their smallest eigenvalue
#'   computed for diagnostics.
#'
#' @return A list with the inverse \code{inv} and the diagnostics
#'   \code{method}, \code{ridge}, \code{min_eig}, \code{n_floored},
#'   \code{size}.
#' @noRd
.ss_invert_block <- function(M, max_tries = 6L, dense_max = 2500L,
                             diag_eigen_max = 1500L) {
  p <- nrow(M)
  if (p == 0L) stop("empty block")
  
  ## Fast path for sparse blocks (2026-10-08): a sparse Cholesky keeps the
  ## inverse as sparse as the relatedness structure allows. Converting a
  ## merged block of unrelated families to dense made every inverse block
  ## full, ~100x slower here and in every later product. Same numbers to
  ## machine precision; anything that is not cleanly PD falls through to the
  ## ridge / eigenvalue-floor path below.
  if (is(M, "sparseMatrix")) {
    Ms <- forceSymmetric((M + t(M)) / 2)
    if (all(is.finite(Ms@x))) {
      out <- tryCatch({
        ch <- Cholesky(Ms, perm = TRUE, LDL = FALSE)
        as(solve(ch, Diagonal(p)), "generalMatrix")
      }, error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(out))
        return(list(inv = out, method = "chol", ridge = 0, min_eig = NA_real_,
                    n_floored = 0L, size = p))
    }
  }
  
  use_dense <- p <= dense_max
  if (use_dense) {
    M <- as.matrix(M)
    M <- (M + t(M)) / 2                       # remove numerical asymmetry
    if (!all(is.finite(M))) stop("block contains NA/NaN/Inf")
    add_ridge <- function(r) M + diag(r, p)
  } else {
    M <- forceSymmetric(M)
    if (!all(is.finite(M@x))) stop("block contains NA/NaN/Inf")
    add_ridge <- function(r) M + Diagonal(p, r)
  }
  
  sc <- mean(abs(diag(M)))
  if (!is.finite(sc) || sc <= 0) sc <- 1
  ridge0 <- 1e-8 * sc
  
  min_eig <- NA_real_
  ridge   <- 0
  for (k in seq_len(max_tries)) {
    Mk  <- if (ridge > 0) add_ridge(ridge) else M
    out <- tryCatch(chol2inv(chol(Mk)), error = function(e) NULL,
                    warning = function(w) NULL)
    if (!is.null(out))
      return(list(inv = out, method = if (ridge > 0) "chol+ridge" else "chol",
                  ridge = ridge, min_eig = min_eig, n_floored = 0L, size = p))
    if (k == 1L && p <= diag_eigen_max)
      min_eig <- tryCatch(min(eigen(as.matrix(M), symmetric = TRUE,
                                    only.values = TRUE)$values),
                          error = function(e) NA_real_)
    ridge <- if (ridge == 0) ridge0 else ridge * 100
  }
  
  e <- eigen(as.matrix(M), symmetric = TRUE)
  lam_max <- max(e$values)
  if (!is.finite(lam_max) || lam_max <= 0)
    stop(sprintf("block unusable: largest eigenvalue = %.4g", lam_max))
  flo <- 1e-8 * lam_max
  nf  <- sum(e$values < flo)
  lam <- pmax(e$values, flo)
  list(inv = e$vectors %*% (t(e$vectors) / lam), method = "eigen-floor",
       ridge = NA_real_, min_eig = min(e$values), n_floored = nf, size = p)
}

#' Invert a matrix block by block in parallel
#'
#' Inverts every block with \code{.ss_invert_block()}, reports failures and
#' regularised blocks by block id, and assembles the block-diagonal inverse in
#' the original index order.
#'
#' @param mat Square matrix to invert.
#' @param merged_blocks List of index vectors partitioning \code{1:nrow(mat)}.
#' @param label Character label used in messages.
#'
#' @return A sparse block-diagonal inverse with a \code{"block_diag"}
#'   attribute summarising how each block was inverted.
#' @noRd
.ss_block_inverse <- function(mat, merged_blocks, label = "matrix") {
  n    <- nrow(mat)
  perm <- unlist(merged_blocks, use.names = FALSE)
  if (length(perm) != n || anyDuplicated(perm) || !all(perm %in% seq_len(n)))
    stop(sprintf("[%s] blocks do not partition 1:%d (%d indices, %d duplicated)",
                 label, n, length(perm), sum(duplicated(perm))))
  
  ## mc.preschedule = FALSE: a worker killed by the OOM reaper then takes down
  ## only its own block, not every block scheduled on that core.
  res <- mclapply(seq_along(merged_blocks), function(i) {
    blk <- merged_blocks[[i]]
    tryCatch(.ss_invert_block(mat[blk, blk, drop = FALSE]),
             error = function(e) structure(
               list(i = i, size = length(blk), msg = conditionMessage(e)),
               class = "ss_block_failure"))
  }, mc.cores = .ss_cores(), mc.preschedule = FALSE)
  
  died   <- which(vapply(res, inherits, logical(1), "try-error"))
  failed <- which(vapply(res, inherits, logical(1), "ss_block_failure"))
  if (length(died) || length(failed)) {
    for (i in died)
      message(sprintf("[%s] block %d (n=%d) WORKER DIED (likely OOM/segfault): %s",
                      label, i, length(merged_blocks[[i]]),
                      conditionMessage(attr(res[[i]], "condition"))))
    for (i in failed)
      message(sprintf("[%s] block %d (n=%d) FAILED: %s",
                      label, i, res[[i]]$size, res[[i]]$msg))
    first <- if (length(failed))
      sprintf("block %d (n=%d): %s", failed[1], res[[failed[1]]]$size, res[[failed[1]]]$msg)
    else sprintf("block %d worker died", died[1])
    stop(sprintf("[%s] %d of %d blocks could not be inverted; first failure -> %s",
                 label, length(died) + length(failed), length(res), first))
  }
  
  meth <- vapply(res, `[[`, character(1), "method")
  bad  <- which(meth != "chol")
  if (length(bad)) {
    detail <- vapply(bad, function(i) sprintf(
      "  [%s] block %d (n=%d): %s | ridge=%.3g | min eigenvalue=%.4g | %d eigenvalue(s) floored",
      label, i, res[[i]]$size, res[[i]]$method, res[[i]]$ridge,
      res[[i]]$min_eig, res[[i]]$n_floored), character(1))
    hdr <- sprintf("[%s] %d of %d blocks were not positive definite",
                   label, length(bad), length(res))
    if (.ss_strict()) stop(paste(c(hdr, detail), collapse = "\n"), call. = FALSE)
    warning(hdr, " and were regularised", call. = FALSE)
    if (.ss_verbose()) for (d in detail) message(d)
  }
  
  out <- bdiag(lapply(res, `[[`, "inv"))
  if (is.unsorted(perm)) {            # bdiag orders by block, not by index
    ip  <- order(perm)
    out <- out[ip, ip, drop = FALSE]
  }
  attr(out, "block_diag") <- data.frame(
    block = seq_along(res),
    size = vapply(res, `[[`, numeric(1), "size"),
    method = meth,
    ridge = vapply(res, function(z) as.numeric(z$ridge), numeric(1)),
    min_eig = vapply(res, function(z) as.numeric(z$min_eig), numeric(1)),
    floored = vapply(res, function(z) as.numeric(z$n_floored), numeric(1)))
  out
}

#' Find Block Structures in a Sparse Matrix
#'
#' This function identifies independent block structures within a sparse matrix
#' (typically a precision or covariance matrix). It uses a vectorized approach
#' to determine the boundaries where the matrix can be partitioned into
#' non-overlapping sub-matrices.
#'
#' Matrices in symmetric storage (e.g. \code{dsCMatrix}) are converted to
#' general storage first, since \code{summary()} on them returns only one
#' triangle.
#'
#' @param A A sparse matrix (usually of class \code{dgCMatrix} or similar).
#' @param threshold Numeric. Elements with absolute values smaller than this
#'   are treated as zero. Default is 0 (use every stored element).
#'
#' @return A list where each element is an integer vector containing the
#'   row/column indices for a specific block.
#' @export
find_blocks_vectorized <- function(A, threshold = 0) {
  n <- nrow(A)
  ## summary() on a symmetric-storage Matrix (dsCMatrix) returns only ONE
  ## triangle. With uplo == "L" every max_j would be <= i, the algorithm would
  ## report n blocks of size 1, and the "block-wise inverse" would silently be
  ## just the inverse of the diagonal. Force general storage first.
  if (is(A, "symmetricMatrix")) A <- as(A, "generalMatrix")
  sm <- summary(A)
  if (threshold > 0) sm <- sm[abs(sm$x) >= threshold, , drop = FALSE]
  
  ## largest non-zero column index in each row (0 if the row is empty)
  max_j <- rep(0L, n)
  if (nrow(sm) > 0) {
    tmp <- tapply(sm$j, sm$i, max)
    max_j[as.integer(names(tmp))] <- as.integer(tmp)
  }
  ## an empty row only affects itself
  zero_rows <- which(max_j == 0)
  if (length(zero_rows) > 0) max_j[zero_rows] <- zero_rows
  
  cum_max <- cummax(max_j)
  block_boundaries <- which(seq_len(n) == cum_max)
  if (!length(block_boundaries) || block_boundaries[length(block_boundaries)] != n)
    block_boundaries <- unique(c(block_boundaries, n))
  
  blocks <- vector("list", length(block_boundaries))
  start <- 1L
  for (i in seq_along(block_boundaries)) {
    end <- block_boundaries[i]
    blocks[[i]] <- start:end
    start <- end + 1L
  }
  blocks
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
  find_blocks_vectorized(A, threshold = threshold)
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
  rbind(cbind(tl, tr), cbind(bl, br))
}

#' Parallel Block-wise Matrix Inversion
#'
#' This function performs an efficient inversion of a large sparse matrix by
#' partitioning it into independent blocks. It merges smaller blocks to
#' optimize parallel processing and uses Cholesky decomposition for
#' numerical stability. A block that is not positive definite is regularised
#' with an escalating ridge, then an eigenvalue-floored inverse, and reported
#' by block id; see \code{\link{SynPALM-package}} for the controlling options.
#'
#' @param wait_matrix A large square matrix (typically a sparse \code{dgCMatrix}) to be inverted.
#' @param threshold Integer. The maximum number of rows/columns for a merged block
#'   to optimize memory and computation. Default is 1000.
#' @param label Character label used in diagnostic messages.
#'
#' @return A sparse block-diagonal matrix representing the inverse of the input
#'   matrix, with a \code{"block_diag"} attribute summarising how each block
#'   was inverted.
#' @export
#' @import Matrix
#' @import parallel
matrix_inv_block <- function(wait_matrix, threshold = 1000, label = "Sigma") {
  blocks        <- find_blocks_vectorized(wait_matrix)
  merged_blocks <- .ss_merge_blocks(blocks, threshold)
  .ss_block_inverse(wait_matrix, merged_blocks, label = label)
}

#' Invert a Relationship Matrix with Adaptive Block Detection
#'
#' This function performs block-wise inversion of a relationship matrix (A-matrix).
#' It features an adaptive thresholding mechanism that increases the sparseness
#' threshold until the maximum block size is manageable (<= \code{max_block}) for inversion.
#'
#' @param Amatrix A large square matrix (typically a sparse genetic relationship matrix).
#' @param thr Numeric. The initial threshold for detecting independent blocks. Default is 0.006.
#' @param threshold Integer. The maximum number of rows/columns for merged blocks
#'   during parallel computation. Default is 1000.
#' @param max_block Integer. \code{thr} is raised until no block is larger than
#'   this. Defaults to \code{threshold}.
#' @param thr_step Numeric. Increment applied to \code{thr} at each step.
#' @param max_thr_iter Integer. Maximum number of increments before giving up
#'   with a warning.
#' @param label Character label used in diagnostic messages.
#'
#' @return A sparse block-diagonal matrix representing the inverse of the
#'   (possibly thresholded) input matrix, with a \code{"block_diag"} attribute
#'   summarising how each block was inverted.
#' @export
#' @import Matrix
#' @import parallel
matrix_inv_Amatrix <- function(Amatrix, thr = 0.006, threshold = 1000,
                               max_block = threshold, thr_step = 0.001,
                               max_thr_iter = 200L, label = "A") {
  it <- 0L
  repeat {
    blocks_A <- find_blocks_vectorized_threshold(Amatrix, threshold = thr)
    if (!length(blocks_A))
      stop(sprintf("[%s] no blocks found at thr = %.4g", label, thr))
    ja <- max(lengths(blocks_A))
    if (ja <= max_block) break
    it <- it + 1L
    if (it > max_thr_iter) {           # the old while() loop had no exit guard
      warning(sprintf("[%s] gave up thresholding at thr = %.4g; largest block still %d",
                      label, thr, ja), call. = FALSE)
      break
    }
    thr <- thr + thr_step
  }
  if (.ss_verbose() && it > 0L)
    message(sprintf("[%s] thr raised to %.4g in %d step(s): %d blocks, largest = %d",
                    label, thr, it, length(blocks_A), ja))
  
  merged_blocks <- .ss_merge_blocks(blocks_A, threshold)
  .ss_block_inverse(Amatrix, merged_blocks, label = label)
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
  score <- as.numeric(Ug^2 / Vu)
  negative_log10_pval <- -pchisq(score, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  list(score = score, negative_log10_pval = negative_log10_pval)
}

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
score_test_SynSurrG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)
  
  snpindex <- 1:ncol(g_matrix)
  
  G_all <- g_matrix #[,snpindex]
  G_obs <- as.matrix(G_all[obs_protein_index,])
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- VV1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(snpindex), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  )
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
score_test_SynSurr_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix <- as.matrix(g_matrix[independent_indices,])
  
  n_obs <- length(obs_protein_index)
  n <- nrow(g_matrix)
  
  snpindex <- 1:ncol(g_matrix)
  
  G_all <- g_matrix #[,snpindex]
  G_obs <- as.matrix(G_all[obs_protein_index,])
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- VV1[snp_d,] %*% G_obs[,snp_d]
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(snpindex), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  )
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
score_test_OracleG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(g_matrix, invS11_res))
  A <- crossprod(SX, g_matrix)
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_OracleG","negative_log10_pval_OracleG")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix))) - as.numeric(colSums(g_matrix * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_Oracle_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix <- as.matrix(g_matrix[independent_indices,])
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(g_matrix, invS11_res))
  A <- crossprod(SX, g_matrix)
  VU <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Oracle","negative_log10_pval_Oracle")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix * (inv_Sigma11 %*% g_matrix))) - as.numeric(colSums(g_matrix * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_ObsG_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  g_matrix_sub <- g_matrix[obs_protein_index,]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) - as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta
  
  return(results)
}

#' Score test of the synthetic phenotype among measured individuals
#'
#' Mixed-model score test for association between each variant and the
#' synthetic phenotype, restricted to individuals who also have a measured
#' target phenotype.
#'
#' @param g_matrix Numeric genotype matrix, individuals in rows and variants
#'   in columns, over the full cohort.
#' @param step1_pars The list returned by
#'   \code{\link{ObsG_hatY_typeIerror_step1}}.
#'
#' @return A named list with the score statistic, \eqn{-\log_{10}} p-value,
#'   effect estimate and its variance, one entry per variant.
#' @export
score_test_ObsG_hatY_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  g_matrix_sub <- g_matrix[obs_protein_index,]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG_hatY","negative_log10_pval_ObsG_hatY")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) - as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_ObsG_hatY <- hat_beta
  results$var_hat_beta_ObsG_hatY <- var_hat_beta
  
  return(results)
}

#' Score test of the synthetic phenotype among unmeasured individuals
#'
#' Mixed-model score test for association between each variant and the
#' synthetic phenotype, restricted to individuals with no measured target
#' phenotype.
#'
#' @param g_matrix Numeric genotype matrix, individuals in rows and variants
#'   in columns, over the full cohort.
#' @param step1_pars The list returned by
#'   \code{\link{unObsG_hatY_typeIerror_step1}}.
#'
#' @return A named list with the score statistic, \eqn{-\log_{10}} p-value,
#'   effect estimate and its variance, one entry per variant.
#' @export
score_test_unObsG_hatY_multiply <- function(g_matrix,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  g_matrix_sub <- g_matrix[unobs_protein_index,]
  
  ## for unobs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_unObsG_hatY","negative_log10_pval_unObsG_hatY")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) - as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_unObsG_hatY <- hat_beta
  results$var_hat_beta_unObsG_hatY <- var_hat_beta
  
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
score_test_Obs_multiply <- function(g_matrix,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix_all_sub <- as.matrix(g_matrix[independent_indices,])
  
  g_matrix_sub <- g_matrix_all_sub[obs_protein_index,]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(g_matrix_sub, invS11_res))
  A <- crossprod(SX, g_matrix_sub)
  VU <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Obs","negative_log10_pval_Obs")
  
  P1G <- X %*% (XtSX_inv %*% A)
  inv_beta <- as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% g_matrix_sub))) - as.numeric(colSums(g_matrix_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_SynSurrG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  n_obs <- length(obs_protein_index)
  n <- length(G_all)
  
  G_obs <- G_all[obs_protein_index]
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- VV1 %*% G_obs
    
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(1), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurrG = surr_result$score,
    negative_log10_pval_SynSurrG = surr_result$negative_log10_pval,
    hat_beta_SynSurrG = hat_beta_syn,
    var_hat_beta_SynSurrG = var_hat_beta_syn
  )
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
score_test_SynSurr_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all <- G_all[independent_indices]
  
  n_obs <- length(obs_protein_index)
  n <- length(G_all)
  
  G_obs <- G_all[obs_protein_index]
  
  temA <- inv_Sigma22 %*% G_all
  A11 <- colSums(G_all * temA) ## t(G) inv_Sigma22 G
  #Att <- inv_Sigma22 %*% X_all
  A12 <- crossprod(G_all, Att) ## t(G) inv_Sigma22 X
  #A22 <- crossprod(X_all, Att)
  
  Stt1 <- crossprod(G_obs, V11) # t(Gobs) V11
  Stt2 <- crossprod(G_obs, Btt1) # t(Gobs) V11 Xobs
  VV1_2 <- Stt2 %*% B_mat # t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV1 <- Stt1 - VV1_2 # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
  
  bt1 <- Sigma12_Sigma22inv %*% G_all
  
  Ah <- t(temA) %*% t(Sigma12)
  
  hh1 <- cbind(rowSums(VV1 * t(bt1)), VV1 %*% bt2)
  
  V1 <- rowSums((VV1 %*% Sigma11) * VV1)
  
  # 定义单个 SNP 计算的函数
  compute_snp <- function(snp_d) {
    
    ## 1. estimate alpha 
    # 构造 A 部分
    A11_sub <- as(as(as(A11[snp_d], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12_sub <- as(as(as(t(A12[snp_d,]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat   <- block_matrix(A11_sub, A12_sub, t(A12_sub), A22) # (Z^T %*% Sigma_22^{-1} %*% Z)
    # 构造 alpha 右端项 
    Ar0 <- rbind(t(temA[,snp_d]),t(Att)) # Z %*% Sigma_22^{-1}
    # 得到 alpha 估计值
    alpha <- solve(A_mat, Ar0 %*% hatY)
    
    
    ## 2. estimate beta
    #Btt1 <- V11 %*% X_obs
    #Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
    #B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    #B1 <- B_mat %*% Y
    
    
    #bt2 <- Sigma12_Sigma22inv %*% X_all
    b1 <- cbind(bt1[,snp_d],bt2)
    #b3 <- Sigma12_Sigma22inv %*% hatY
    
    #Sigma12_Sigma22inv_rs_alpha <- b3 - b1 %*% alpha # Sigma12_Sigma22inv (hatY-Z alpha)
    #beta <- B1 - B_mat %*% Sigma12_Sigma22inv_rs_alpha
    beta <- B1 - B2 + B_mat %*% (b1 %*% alpha)
    
    ## 3. calculate the score
    S_syn <- as.numeric(Stt1[snp_d,] %*% (Y - X_obs %*% beta - B2_1 + bt1[,snp_d] * alpha[1] + bt2 %*% alpha[-1]))
    
    ## 4. calculate the variance of the score
    
    #the first element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV1 <- t(Stt1[snp_d,] - VV1_2[snp_d,]) # t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
    
    
    #the second element of (t(Gobs) V11 - t(Gobs) V11 Xobs (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11) %*% M de (first part)
    #VV2
    #VV2_1 <- -VV1 %*% Sigma12_Sigma22inv %*% t(Sigma12)
    Arr <- rbind(Ah[snp_d,],Atb)
    Ar <- solve(A_mat, Arr)
    
    VV2 <- VV2_1[snp_d,] + hh1[snp_d,] %*% Ar
    
    V2 <- VV2 %*% VV1[snp_d,]
    
    #V1 <- VV1 %*% Sigma11  %*% t(VV1)
    #V2 <- VV2 %*% t(Sigma12)  %*% t(VV1)
    #V3 <- VV2 %*% Sigma22  %*% t(VV2)
    
    VU <- as.numeric(V1[snp_d] + V2)
    
    inv_beta <- VV1 %*% G_obs
    hat_beta <- as.numeric(S_syn/inv_beta)
    var_hat_beta <- as.numeric(VU/(inv_beta)^2)
    
    list(#alpha = alpha, beta = beta, 
      S_syn=S_syn, VU = VU, hat_beta = hat_beta, var_hat_beta = var_hat_beta)
  }
  
  # 利用 lapply 对所有 SNP 进行计算
  results <- lapply(seq_along(1), compute_snp)
  
  # 将各个结果整合成矩阵和向量
  #hat_alpha   <- do.call(cbind, lapply(results, `[[`, "alpha"))
  #hat_beta    <- do.call(cbind, lapply(results, `[[`, "beta"))
  VU_syn <- sapply(results, `[[`, "VU")
  S_syn <- sapply(results, `[[`, "S_syn")
  hat_beta_syn <- sapply(results, `[[`, "hat_beta")
  var_hat_beta_syn <- sapply(results, `[[`, "var_hat_beta")
  
  #S_obs <- as.numeric(crossprod(G_obs, invS11_res))
  #A <- crossprod(SX, G_obs)
  #VU_obs <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * XtSX_inv %*% A))
  
  surr_result <- compute_score(S_syn, VU_syn)
  #obs_result <- compute_score(S_obs, VU_obs)
  
  list(
    T_score_SynSurr = surr_result$score,
    negative_log10_pval_SynSurr = surr_result$negative_log10_pval,
    hat_beta_SynSurr = hat_beta_syn,
    var_hat_beta_SynSurr = var_hat_beta_syn
  )
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
score_test_OracleG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(G_all, invS11_res))
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_OracleG","negative_log10_pval_OracleG")
  
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all)
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all))) - as.numeric(colSums(G_all * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_Oracle_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all <- G_all[independent_indices]
  
  ## for oracle test statistic
  SU <- as.numeric(crossprod(G_all, invS11_res))
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Oracle","negative_log10_pval_Oracle")
  
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all)
  inv_beta <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all))) - as.numeric(colSums(G_all * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_ObsG_single <- function(G_all,step1_pars) {
  
  list2env(step1_pars, envir = environment())
  G_all_sub <- G_all[obs_protein_index]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(G_all_sub, invS11_res))
  A <- crossprod(SX, G_all_sub)
  VU <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_all_sub)
  inv_beta <- as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% G_all_sub))) - as.numeric(colSums(G_all_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
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
score_test_Obs_single <- function(G_all,step1_pars,independent_indices) {
  
  list2env(step1_pars, envir = environment())
  
  G_all_sub <- G_all[independent_indices]
  G_sub <- G_all_sub[obs_protein_index]
  
  ## for obs test statistic
  SU <- as.numeric(crossprod(G_sub, invS11_res))
  A <- crossprod(SX, G_sub)
  VU <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_Obs","negative_log10_pval_Obs")
  
  P1G <- X %*% (XtSX_inv %*% t(SX) %*% G_sub)
  inv_beta <- as.numeric(colSums(G_sub * (inv_Sigma11 %*% G_sub))) - as.numeric(colSums(G_sub * (inv_Sigma11 %*% P1G)))
  hat_beta <- SU/(inv_beta)
  var_hat_beta <- VU/(inv_beta)^2
  
  results$hat_beta_Obs <- hat_beta
  results$var_hat_beta_Obs <- var_hat_beta
  
  return(results)
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
DGP <- function(n_obs,miss,rho,maf = 0.25,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5, 
                pve_g = 0, pve_x = 0.10) {
  
  ## using rho to generate the tauTS and sigmaTS
  tauTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*3/5
  sigmaTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*2/5
  
  # Subjects with observed target outcomes.
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype and covariate.
  G_all <- stats::rbinom(n = n_all, size = 2, prob = maf)
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all, X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - pve_g - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- sqrt(pve_g * V_total /var(G_all))
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # generate the block-wise GRM matrix
  blocks <- seq(1, n_all-1, by = 2)
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)  # 行索引
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)  # 列索引
  values <- rep(c(1, 0.5, 0.5, 1), times = length(blocks))  # 元素值
  
  GRM <- sparseMatrix(i = rows,j = cols,x = values, dims = c(n_all, n_all))
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  ## generate observed Y_obs and suggrate S
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)
  
  I_matrix <- sparseMatrix(
    i = 1:n_obs,j = 1:n_obs,x = sigmaTS,
    dims = c(n_obs, n_all)
  )
  
  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)  # 对称
  
  # Sigma22 的逆
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  mu_cond <- mu_all[1:n_obs] + Sigma12 %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond <- Sigma11 - Sigma12 %*% inv_Sigma22 %*% Sigma21
  #L_cond <- chol(Sigma_cond)
  #Sigma_cond <- as.matrix(nearPD(Sigma_cond)$mat)
  L_cond <- chol(Sigma_cond)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )
  
  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)  # 对称
  
  Sigma12_oracle_Sigma22inv <- Sigma12_oracle %*% inv_Sigma22
  
  mu_cond_oracle <- mu_all + Sigma12_oracle %*% inv_Sigma22 %*% (S - mu_all)
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle %*% inv_Sigma22 %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)
  
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
DGP_GRMr <- function(n_obs,miss,GRMr){
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # generate the block-wise GRM matrix
  blocks <- seq(1, n_all-1, by = 2)
  rows <- rep(blocks, each = 4) + c(0, 0, 1, 1)  # 行索引
  cols <- rep(blocks, each = 4) + c(0, 1, 0, 1)  # 列索引
  values <- rep(c(1, GRMr, GRMr, 1), times = length(blocks))  # 元素值
  
  GRM <- sparseMatrix(i = rows,j = cols,x = values, dims = c(n_all, n_all))
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
DGP_GRMr_nullmodel_step1 <- function(n_obs,miss,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,GRM){
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  ## 预先计算一些量，避免重复
  Sigma11 <- tauT^2 * GRM_obs + sigmaT^2 * Diagonal(n = n_obs, x = 1)
  Sigma22 <- tauS^2 * GRM + sigmaS^2 * Diagonal(n = n_all, x = 1)
  
  # Sigma22 的逆
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  Sigma11_oracle <- tauT^2 * GRM + sigmaT^2 * Diagonal(n = n_all, x = 1)
  
  params <- c("Sigma11","Sigma22","inv_Sigma22","Sigma11_oracle","chol_Sigma22")
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
DGP_GRMr_nullmodel_step2 <- function(n_obs,miss,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,
                                     rho,GRM,pre_vars_step1){
  
  list2env(pre_vars_step1, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  GRM_obs <- GRM[1:n_obs,1:n_obs]
  GRM_obs_all <- GRM[1:n_obs,]
  colnames(GRM_obs) <- 1:n_obs
  rownames(GRM_obs) <- 1:n_obs
  
  colnames(GRM) <- 1:n_all
  rownames(GRM) <- 1:n_all
  
  ## using rho to generate the tauTS and sigmaTS
  tauTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*3/5
  sigmaTS <- sqrt((tauT^2+sigmaT^2)*(tauS^2+sigmaS^2))*rho*2/5
  
  I_matrix <- sparseMatrix(
    i = 1:n_obs,j = 1:n_obs,x = sigmaTS,
    dims = c(n_obs, n_all)
  )
  
  Sigma12 <- tauTS * GRM_obs_all + I_matrix
  Sigma21 <- t(Sigma12)  # 对称
  
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  
  Sigma_cond <- Sigma11 - Sigma12_Sigma22inv %*% Sigma21
  #L_cond <- chol(Sigma_cond)
  #Sigma_cond <- as.matrix(nearPD(Sigma_cond)$mat)
  L_cond <- chol(Sigma_cond)
  
  I_matrix_oracle <- sparseMatrix(
    i = 1:n_all, j = 1:n_all, x = sigmaTS, dims = c(n_all, n_all)
  )
  
  Sigma12_oracle <- tauTS * GRM + I_matrix_oracle
  Sigma21_oracle <- t(Sigma12_oracle)  # 对称
  
  Sigma12_oracle_Sigma22inv <- Sigma12_oracle %*% inv_Sigma22
  
  Sigma_cond_oracle <- Sigma11_oracle - Sigma12_oracle_Sigma22inv %*% Sigma21_oracle
  L_cond_oracle <- chol(Sigma_cond_oracle)
  
  params <- c("chol_Sigma22","Sigma12_Sigma22inv","L_cond",
              "Sigma12_oracle_Sigma22inv","L_cond_oracle")
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
DGP_GRMr_nullmodel_step3 <- function(n_obs,miss,tauT = 0.7, sigmaT = 0.7, pve_x = 0.10, 
                                     pre_vars_step2, GRM) {
  
  list2env(pre_vars_step2, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype 0 and covariate.
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(rep(0,n_all),X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - 0 - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- 0
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  # 生成hatY
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  # 生成Y
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)
  
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
typeIerror_generate_gmatrix <- function(n_obs,miss,chunk_size,maf){
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  g_matrix <- matrix(rbinom(n_all * chunk_size, size = 2, prob = maf), nrow = n_all, ncol = chunk_size)
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
generate_gmatrix_kinship <- function(n_obs,miss,maf,kinsetting){
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  draw_gamete <- function(g) rbinom(length(g), size = 1, prob = g / 2)
  
  g_grand_father <- rbinom(n_all/2, size = 2, prob = maf)
  g_grand_mother <- rbinom(n_all/2, size = 2, prob = maf)
  
  if(kinsetting == "unkin"){
    geno_vec <- c(rbind(g_grand_father, g_grand_mother))
  }
  if(kinsetting == "siblings"){
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    geno_vec <- c(rbind(g_father, g_uncle))
  }
  if(kinsetting == "cousin"){
    g_father <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_uncle <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_mother <- rbinom(n_all/2, size = 2, prob = maf)
    g_me <- draw_gamete(g_father) + draw_gamete(g_mother)
    g_rand <- rbinom(n_all/2, size = 2, prob = maf)
    g_far_sister <- draw_gamete(g_uncle) + draw_gamete(g_rand)
    geno_vec <- c(rbind(g_me, g_far_sister))
  }
  if(kinsetting == "halfsiblings"){
    g_uncle1 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother)
    g_grand_mother2 <- rbinom(n_all/2, size = 2, prob = maf)
    g_uncle2 <- draw_gamete(g_grand_father) + draw_gamete(g_grand_mother2)
    geno_vec <- c(rbind(g_uncle1, g_uncle2))
  }
  return(geno_vec)
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
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  # r <- rank(test_df %>% pull(Y_all))
  # test_df$Y_all_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of y_all
  # test_df <- INT(test_df, "Y_obs") 
  # r <- rank(test_df %>% pull(yhat))
  # test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat
  # 
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )%>% na.omit()
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )
  
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurrG)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurrG, na.action = na.exclude)
  
  ## residuals aligned to their index sets. na.action = na.exclude keeps the NA
  ## padding, so res_S[obs_protein_index] refers to the right individuals even
  ## if lm() dropped rows; the old code indexed a shortened residual vector by
  ## position, which silently misaligns whenever a covariate has an NA.
  res_T <- as.numeric(model_lm1$residuals)
  res_S <- as.numeric(residuals(model_lm2))
  if (length(res_T) != n_obs)
    stop("SynSurrG_typeIerror_step1: lm1 residual length ", length(res_T), " != n_obs ", n_obs,
         " -- NA in covariates or surrogate at observed individuals")
  if (length(res_S) != n)
    stop("SynSurrG_typeIerror_step1: lm2 residual length ", length(res_S), " != n ", n)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(res_T) %*% GRM_o %*% res_T)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_typeIerror_step1")
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  ## inv_Sigma11 is never used downstream (it is not in `params`): dropping it
  ## saves one full block-wise inversion per protein and removes a crash site.
  ## inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat", 
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12",
    "vc"
  )
  
  SynSurrG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(SynSurrG_typeIerror_step1_pars)
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
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )%>% na.omit()
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )
  
  preds <- grep("^x\\.X\\.", names(reml_data_SynSurrG), value = TRUE)
  fml <- reformulate(preds, response = "y")
  
  model_lm1 <- lm(fml, data = reml_data_SynSurrG)
  
  preds <- grep("^x\\.X\\.", names(lm_data_SynSurrG), value = TRUE)
  fml <- reformulate(preds, response = "haty")
  model_lm2 <- lm(fml, data = lm_data_SynSurrG, na.action = na.exclude)
  
  ## residuals aligned to their index sets. na.action = na.exclude keeps the NA
  ## padding, so res_S[obs_protein_index] refers to the right individuals even
  ## if lm() dropped rows; the old code indexed a shortened residual vector by
  ## position, which silently misaligns whenever a covariate has an NA.
  res_T <- as.numeric(model_lm1$residuals)
  res_S <- as.numeric(residuals(model_lm2))
  if (length(res_T) != n_obs)
    stop("SynSurrG_ablation_estimate: lm1 residual length ", length(res_T), " != n_obs ", n_obs,
         " -- NA in covariates or surrogate at observed individuals")
  if (length(res_S) != n)
    stop("SynSurrG_ablation_estimate: lm2 residual length ", length(res_S), " != n ", n)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(res_T) %*% GRM_o %*% res_T)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_ablation_estimate")
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  ## inv_Sigma11 is never used downstream (it is not in `params`): dropping it
  ## saves one full block-wise inversion per protein and removes a crash site.
  ## inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df[, grepl("^X\\.", names(test_df)), drop = FALSE])
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12",
    "vc"
  )
  
  SynSurrG_ablation_est_pars <- mget(params, envir = environment())
  
  return(SynSurrG_ablation_est_pars)
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
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  # r <- rank(test_df %>% pull(Y_all))
  # test_df$Y_all_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of y_all
  # test_df <- INT(test_df, "Y_obs") 
  # r <- rank(test_df %>% pull(yhat))
  # test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat
  
  reml_data_oracle <- data.frame(
    y = test_df$Y_all,
    x = test_df$X
  )
  
  model_lm <- lm(y ~ x, data = reml_data_oracle)
  
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "OracleG_typeIerror_step1")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)
  
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- reml_data_oracle$y
  X <- cbind(rep(1,n),reml_data_oracle$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
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
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  reml_data_oracle <- data.frame(
    y = test_df$Y_all,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )
  
  preds <- grep("^x\\.X\\.", names(reml_data_oracle), value = TRUE)
  fml <- reformulate(preds, response = "y")
  
  model_lm <- lm(fml, data = reml_data_oracle)
  
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "OracleG_ablation_estimate")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)
  
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- reml_data_oracle$y
  X <- as.matrix(cbind(rep(1,n),test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]))
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  OracleG_ablation_est_pars <- mget(params, envir = environment())
  
  return(OracleG_ablation_est_pars)
  
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
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  
  data_obs <- data.frame(
    y = test_df$Y_obs,
    x = test_df$X
  )%>% na.omit()
  
  n_obs <- nrow(data_obs)
  
  model_lm <- lm(y ~ x, data = data_obs)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_typeIerror_step1")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- cbind(rep(1,length(Y)),data_obs$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  ObsG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(ObsG_typeIerror_step1_pars)
  
}

#' Null model for the synthetic phenotype among measured individuals
#'
#' Fits the linear mixed null model to the synthetic (predicted) phenotype,
#' restricted to individuals who also have a measured target phenotype, and
#' returns the quantities the corresponding score test reuses across variants.
#' Used to check that the synthetic phenotype alone is well calibrated.
#'
#' @param mydf A list with elements \code{X_all} (covariates for all
#'   individuals, columns prefixed "X."), \code{Y} (oracle target
#'   phenotype), \code{Y_obs} (measured target phenotype, \code{NA} where
#'   unmeasured), \code{S} (synthetic phenotype) and \code{GRM} (sparse
#'   genetic relatedness matrix).
#'
#' @return A named list of precomputed step-1 quantities, to be passed as
#'   \code{step1_pars} to \code{\link{score_test_ObsG_hatY_multiply}}.
#' @export
ObsG_hatY_typeIerror_step1 <- function(mydf) {
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  
  data_obs <- data.frame(
    y = test_df$yhat[obs_protein_index],
    x = test_df$X[obs_protein_index])
  
  n_obs <- nrow(data_obs)
  
  model_lm <- lm(y ~ x, data = data_obs)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_hatY_typeIerror_step1")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- cbind(rep(1,length(Y)),data_obs$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  ObsG_hatY_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(ObsG_hatY_typeIerror_step1_pars)
  
}

#' Null model for the synthetic phenotype among unmeasured individuals
#'
#' As \code{\link{ObsG_hatY_typeIerror_step1}}, but restricted to the
#' individuals with no measured target phenotype. Fitting the two subsets
#' separately makes it possible to check that the synthetic phenotype behaves
#' the same way in both.
#'
#' @param mydf A list with elements \code{X_all}, \code{Y}, \code{Y_obs},
#'   \code{S} and \code{GRM}; see \code{\link{ObsG_hatY_typeIerror_step1}}.
#'
#' @return A named list of precomputed step-1 quantities, to be passed as
#'   \code{step1_pars} to \code{\link{score_test_unObsG_hatY_multiply}}.
#' @export
unObsG_hatY_typeIerror_step1 <- function(mydf) {
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  unobs_protein_index <- which(is.na(test_df$Y_obs))
  GRM_unobs_unobs <- GRM[unobs_protein_index, unobs_protein_index]
  
  
  data_unobs <- data.frame(
    y = test_df$yhat[unobs_protein_index],
    x = test_df$X[unobs_protein_index])
  
  n_unobs <- nrow(data_unobs)
  
  model_lm <- lm(y ~ x, data = data_unobs)
  
  GRM_uno <- GRM_unobs_unobs
  diag(GRM_uno) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_uno %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_uno)^2) #equal to sum(diag((GRM_uno) %*% t(GRM_unobs_unobs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_unobs_unobs)))/nrow(GRM_unobs_unobs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "unObsG_hatY_typeIerror_step1")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM_unobs_unobs + hat_sigma_T2 * Diagonal(n_unobs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_unobs$y
  X <- cbind(rep(1,length(Y)),data_unobs$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "unobs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  unObsG_hatY_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(unObsG_hatY_typeIerror_step1_pars)
  
}

#' Ordinary least squares null model for the measured target phenotype
#'
#' Fits the null model by ordinary least squares, ignoring relatedness, to the
#' individuals with a measured target phenotype. Provides the reference
#' against which the mixed-model analyses are compared: any inflation this
#' shows but the mixed model does not is attributable to relatedness.
#'
#' @param mydf A list with elements \code{X_all}, \code{Y}, \code{Y_obs}
#'   and \code{S}. The genetic relatedness matrix is not used.
#'
#' @return A named list holding the residuals, residual variance estimate and
#'   \eqn{(X'X)^{-1}}, to be passed as \code{step1_pars} to
#'   \code{\link{score_test_ObsG_lm_multiply}}.
#' @export
ObsG_typeIerror_step1_lm <- function(mydf) {
  
  test_df <- data.frame(
    sample_index = seq_along(mydf$Y_obs),
    X = mydf$X_all,
    Y_all = mydf$Y,
    Y_obs = mydf$Y_obs,
    yhat = mydf$S
  )
  
  ## observed protein samples
  data_obs <- test_df[complete.cases(test_df[, c("Y_obs", "X")]), ]
  
  obs_protein_index <- data_obs$sample_index
  
  Y <- data_obs$Y_obs
  
  ## X needs to include intercept
  X <- cbind(
    Intercept = 1,
    x = data_obs$X
  )
  
  n_obs <- length(Y)
  p <- ncol(X)
  
  ## ordinary least squares under the null model: Y ~ X
  XtX_inv <- solve(crossprod(X))
  alpha_hat <- XtX_inv %*% crossprod(X, Y)
  
  residual <- as.numeric(Y - X %*% alpha_hat)
  
  ## residual variance estimate
  hat_sigma2 <- sum(residual^2) / (n_obs - p)
  
  ## quantities used in ordinary LM score test
  XtX_inv <- solve(crossprod(X))
  
  params <- c(
    "obs_protein_index",
    "Y",
    "X",
    "residual",
    "hat_sigma2",
    "XtX_inv",
    "n_obs"
  )
  
  ObsG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(ObsG_typeIerror_step1_pars)
}

#' Ordinary least squares null model for the synthetic phenotype
#'
#' As \code{\link{ObsG_typeIerror_step1_lm}}, but fitted to the synthetic
#' phenotype among the individuals with no measured target phenotype.
#'
#' @param mydf A list with elements \code{X_all}, \code{Y}, \code{Y_obs}
#'   and \code{S}. The genetic relatedness matrix is not used.
#'
#' @return A named list holding the residuals, residual variance estimate and
#'   \eqn{(X'X)^{-1}}, to be passed as \code{step1_pars} to
#'   \code{\link{score_test_unObsG_hatY_lm_multiply}}.
#' @export
unObsG_hatY_typeIerror_step1_lm <- function(mydf) {
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  unobs_protein_index <- which(is.na(test_df$Y_obs))
  
  data_unobs <- data.frame(
    y = test_df$yhat[unobs_protein_index],
    x = test_df$X[unobs_protein_index])
  
  
  Y <- data_unobs$y
  
  ## X needs to include intercept
  X <- cbind(
    Intercept = 1,
    x = data_unobs$x
  )
  
  n_unobs <- length(Y)
  p <- ncol(X)
  
  ## ordinary least squares under the null model: Y ~ X
  XtX_inv <- solve(crossprod(X))
  alpha_hat <- XtX_inv %*% crossprod(X, Y)
  
  residual <- as.numeric(Y - X %*% alpha_hat)
  
  ## residual variance estimate
  hat_sigma2 <- sum(residual^2) / (n_unobs - p)
  
  ## quantities used in ordinary LM score test
  XtX_inv <- solve(crossprod(X))
  
  params <- c(
    "unobs_protein_index",
    "Y",
    "X",
    "residual",
    "hat_sigma2",
    "XtX_inv",
    "n_unobs"
  )
  
  unObsG_hatY_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(unObsG_hatY_typeIerror_step1_pars)
}

#' Ordinary least squares score test, measured individuals
#'
#' Score test under ordinary linear regression, ignoring relatedness, for the
#' measured target phenotype.
#'
#' @param g_matrix Numeric genotype matrix, individuals in rows and variants
#'   in columns, over the full cohort.
#' @param step1_pars The list returned by
#'   \code{\link{ObsG_typeIerror_step1_lm}}.
#'
#' @return A named list with the score statistic, \eqn{-\log_{10}} p-value,
#'   effect estimate and its variance, one entry per variant.
#' @export
score_test_ObsG_lm_multiply <- function(g_matrix, step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix_sub <- as.matrix(g_matrix[obs_protein_index, , drop = FALSE])
  
  ## score numerator: G^T residual
  S <- as.numeric(crossprod(g_matrix_sub, residual))
  
  ## residualize G with respect to X
  A <- crossprod(X, g_matrix_sub)
  
  G_MX_G <- as.numeric(
    colSums(g_matrix_sub * g_matrix_sub) -
      colSums(A * (XtX_inv %*% A))
  )
  
  ## score test under ordinary linear regression
  SU <- S / hat_sigma2
  VU <- G_MX_G / hat_sigma2
  
  results <- compute_score(SU, VU)
  
  names(results) <- c(
    "T_score_ObsG",
    "negative_log10_pval_ObsG"
  )
  
  ## beta estimate and variance
  hat_beta <- S / G_MX_G
  var_hat_beta <- hat_sigma2 / G_MX_G
  
  results$hat_beta_ObsG <- hat_beta
  results$var_hat_beta_ObsG <- var_hat_beta
  
  return(results)
}

#' Ordinary least squares score test, synthetic phenotype
#'
#' Score test under ordinary linear regression, ignoring relatedness, for the
#' synthetic phenotype among individuals with no measured target phenotype.
#'
#' @param g_matrix Numeric genotype matrix, individuals in rows and variants
#'   in columns, over the full cohort.
#' @param step1_pars The list returned by
#'   \code{\link{unObsG_hatY_typeIerror_step1_lm}}.
#'
#' @return A named list with the score statistic, \eqn{-\log_{10}} p-value,
#'   effect estimate and its variance, one entry per variant.
#' @export
score_test_unObsG_hatY_lm_multiply <- function(g_matrix, step1_pars) {
  
  list2env(step1_pars, envir = environment())
  
  g_matrix_sub <- as.matrix(g_matrix[unobs_protein_index, , drop = FALSE])
  
  ## score numerator: G^T residual
  S <- as.numeric(crossprod(g_matrix_sub, residual))
  
  ## residualize G with respect to X
  A <- crossprod(X, g_matrix_sub)
  
  G_MX_G <- as.numeric(
    colSums(g_matrix_sub * g_matrix_sub) -
      colSums(A * (XtX_inv %*% A))
  )
  
  ## score test under ordinary linear regression
  SU <- S / hat_sigma2
  VU <- G_MX_G / hat_sigma2
  
  results <- compute_score(SU, VU)
  
  names(results) <- c(
    "T_score_unObsG_hatY",
    "negative_log10_pval_ObsG_hatY"
  )
  
  ## beta estimate and variance
  hat_beta <- S / G_MX_G
  var_hat_beta <- hat_sigma2 / G_MX_G
  
  results$hat_beta_unObsG_hatY <- hat_beta
  results$var_hat_beta_unObsG_hatY <- var_hat_beta
  
  return(results)
}

#' SynPALM null model with the tested variant retained in the surrogate model
#'
#' As \code{\link{SynSurrG_ablation_estimate}}, except that the variant under
#' test is kept as a covariate when the synthetic phenotype is regressed on the
#' covariates. This is the appropriate step 1 when the prediction model that
#' generated the synthetic phenotype already contains that variant, so that the
#' variance components are not contaminated by the signal being tested. Because
#' step 1 then depends on the variant, it must be refitted for each one, which
#' is far more costly than the shared-step-1 route.
#'
#' @param mydf A list with elements \code{X_all}, \code{Y}, \code{Y_obs},
#'   \code{S} and \code{GRM}.
#' @param G Numeric genotype vector for the single variant to be retained as a
#'   covariate in the surrogate model.
#'
#' @return A named list of precomputed step-1 quantities, to be passed as
#'   \code{step1_pars} to \code{\link{score_test_SynSurrG_single}}.
#' @export
SynSurrG_eachG_step1 <- function(mydf,G) {
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S, G=G)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  # r <- rank(test_df %>% pull(Y_all))
  # test_df$Y_all_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of y_all
  # test_df <- INT(test_df, "Y_obs") 
  # r <- rank(test_df %>% pull(yhat))
  # test_df$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat
  # 
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )%>% na.omit()
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE],
    G = test_df$G
  )
  
  preds <- grep("^x\\.X\\.", names(reml_data_SynSurrG), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm1 <- lm(fml, data = reml_data_SynSurrG)
  
  preds <- c(grep("^x\\.X\\.", names(lm_data_SynSurrG), value = TRUE),"G")
  fml <- reformulate(preds, response = "haty")
  model_lm2 <- lm(fml, data = lm_data_SynSurrG, na.action = na.exclude)
  
  ## residuals aligned to their index sets. na.action = na.exclude keeps the NA
  ## padding, so res_S[obs_protein_index] refers to the right individuals even
  ## if lm() dropped rows; the old code indexed a shortened residual vector by
  ## position, which silently misaligns whenever a covariate has an NA.
  res_T <- as.numeric(model_lm1$residuals)
  res_S <- as.numeric(residuals(model_lm2))
  if (length(res_T) != n_obs)
    stop("SynSurrG_eachG_step1: lm1 residual length ", length(res_T), " != n_obs ", n_obs,
         " -- NA in covariates or surrogate at observed individuals")
  if (length(res_S) != n)
    stop("SynSurrG_eachG_step1: lm2 residual length ", length(res_S), " != n ", n)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(res_T) %*% GRM_o %*% res_T)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_eachG_step1")
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  ## inv_Sigma11 is never used downstream (it is not in `params`): dropping it
  ## saves one full block-wise inversion per protein and removes a crash site.
  ## inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df[, grepl("^X\\.", names(test_df)), drop = FALSE])
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12",
    "vc"
  )
  
  
  SynSurrG_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(SynSurrG_typeIerror_step1_pars)
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
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  data_obs <- data.frame(
    y = test_df$Y_obs,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )%>% na.omit()
  
  n_obs <- nrow(data_obs)
  
  preds <- grep("^x\\.X\\.", names(data_obs), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm <- lm(fml, data = data_obs)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_ablation_estimate")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- as.matrix(cbind(rep(1,length(Y)),data_obs[, grepl("^x\\.X\\.", names(data_obs)), drop = FALSE]))
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  
  
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
SynSurr_typeIerror_step1 <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  )%>% na.omit()
  
  lm_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent$X
  )
  
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurr)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurr)
  
  hat_sigma_T2 <- var(model_lm1$residuals)
  hat_sigma_S2 <- var(model_lm2$residuals)
  hat_sigma_TS <- cov(model_lm2$residuals[obs_protein_index],model_lm1$residuals)
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  chol_A <- chol(A)         # Cholesky
  V11 <- chol2inv(chol_A) 
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df_independent$Y_obs[obs_protein_index]
  hatY <- test_df_independent$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df_independent$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12"
  )
  
  SynSurr_typeIerror_step1_pars <- mget(params, envir = environment())
  
  return(SynSurr_typeIerror_step1_pars)
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
Oracle_typeIerror_step1 <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  data_oracle <- data.frame(
    y = test_df_independent$Y_all,
    x = test_df_independent$X
  )
  
  model_lm <- lm(y ~ x, data = data_oracle)
  
  hat_sigma_T2 <- var(model_lm$residuals)
  
  Sigma11 <- hat_sigma_T2 * Diagonal(n)
  
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_oracle$y
  X <- cbind(rep(1,n),data_oracle$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv","X"
  )
  
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
Obs_typeIerror_step1 <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  data_obs <- data.frame(
    y = test_df_independent$Y_obs,
    x = test_df_independent$X
  )%>% na.omit()
  
  n_obs <- nrow(data_obs)
  
  model_lm <- lm(y ~ x, data = data_obs)
  
  hat_sigma_T2 <- var(model_lm$residuals)
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- cbind(rep(1,length(Y)),data_obs$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X"
  )
  
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
SynSurr_ablation_estimate <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices,], 
                                    Y_all = mydf$Y[independent_indices], Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE]
  )%>% na.omit()
  
  lm_data_SynSurr <- data.frame(
    y = test_df_independent$Y_obs,
    haty = test_df_independent$yhat,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE]
  )
  
  preds <- grep("^x\\.X\\.", names(reml_data_SynSurr), value = TRUE)
  fml <- reformulate(preds, response = "y")
  model_lm1 <- lm(fml, data = reml_data_SynSurr)
  
  preds <- grep("^x\\.X\\.", names(lm_data_SynSurr), value = TRUE)
  fml <- reformulate(preds, response = "haty")
  model_lm2 <- lm(fml, data = lm_data_SynSurr)
  
  hat_sigma_T2 <- var(model_lm1$residuals)
  hat_sigma_S2 <- var(model_lm2$residuals)
  hat_sigma_TS <- cov(model_lm2$residuals[obs_protein_index],model_lm1$residuals)
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  chol_Sigma22 <- chol(Sigma22)         # Cholesky
  inv_Sigma22 <- chol2inv(chol_Sigma22) 
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  chol_A <- chol(A)         # Cholesky
  V11 <- chol2inv(chol_A) 
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df_independent$Y_obs[obs_protein_index]
  hatY <- test_df_independent$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE])
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12"
  )
  
  SynSurr_ablation_est_pars <- mget(params, envir = environment())
  
  return(SynSurr_ablation_est_pars)
}

#' Oracle null model on independent samples
#'
#' Step 1 for the oracle analysis, which uses the complete target phenotype
#' that would be available if nobody were missing a measurement. Restricted to
#' one individual per relatedness block, so the covariance is diagonal and no
#' genetic relatedness matrix is needed. Provides the upper bound on power in
#' the ablation comparisons.
#'
#' @param mydf A list with elements \code{X_all} (covariates, columns prefixed
#'   "X."), \code{Y} (complete target phenotype), \code{Y_obs} and \code{S}.
#' @param independent_indices Integer vector of row indices giving one
#'   individual per relatedness block, as returned by sampling within the
#'   output of \code{\link{find_blocks_vectorized}}.
#'
#' @return A named list of precomputed step-1 quantities, to be passed as
#'   \code{step1_pars} to \code{\link{score_test_Oracle_multiply}}.
#' @export
Oracle_ablation_estimate <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices,], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  n <- nrow(test_df_independent)
  
  data_oracle <- data.frame(
    y = test_df_independent$Y_all,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE]
  )
  
  preds <- grep("^x\\.X\\.", names(data_oracle), value = TRUE)
  fml <- reformulate(preds, response = "y")
  
  model_lm <- lm(fml, data = data_oracle)
  
  hat_sigma_T2 <- var(model_lm$residuals)
  
  Sigma11 <- hat_sigma_T2 * Diagonal(n)
  
  chol_Sigma11 <- chol(Sigma11)         # Cholesky
  inv_Sigma11 <- chol2inv(chol_Sigma11) 
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_oracle$y
  X <- as.matrix(cbind(rep(1,n),test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE]))
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX) 
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "invS11_res","SX","inv_Sigma11","XtSX_inv","X"
  )
  
  Oracle_ablation_est_pars <- mget(params, envir = environment())
  
  return(Oracle_ablation_est_pars)
  
}

#' Observed-only null model on independent samples
#'
#' Step 1 for the observed-only analysis, which discards the synthetic
#' phenotype and uses just the individuals with a measured target phenotype.
#' Restricted to one individual per relatedness block, so the covariance is
#' diagonal and no genetic relatedness matrix is needed. Provides the baseline
#' against which the power gain from the synthetic phenotype is measured.
#'
#' @param mydf A list with elements \code{X_all} (covariates, columns prefixed
#'   "X."), \code{Y}, \code{Y_obs} (\code{NA} where unmeasured) and \code{S}.
#' @param independent_indices Integer vector of row indices giving one
#'   individual per relatedness block, as returned by sampling within the
#'   output of \code{\link{find_blocks_vectorized}}.
#'
#' @return A named list of precomputed step-1 quantities, to be passed as
#'   \code{step1_pars} to \code{\link{score_test_Obs_multiply}}.
#' @export
Obs_ablation_estimate <- function(mydf,independent_indices) {
  
  test_df_independent <- data.frame(X = mydf$X_all[independent_indices,], 
                                    Y_all = mydf$Y[independent_indices], 
                                    Y_obs = mydf$Y_obs[independent_indices], 
                                    yhat = mydf$S[independent_indices])
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df_independent$Y_obs))
  
  data_obs <- data.frame(
    y = test_df_independent$Y_obs,
    x = test_df_independent[, grepl("^X\\.", names(test_df_independent)), drop = FALSE]
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
  X <- as.matrix(cbind(rep(1,length(Y)),data_obs[, grepl("^x\\.X\\.", names(data_obs)), drop = FALSE]))
  
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
DGP_GRMr_gmodel_step3 <- function(n_obs,miss,tauT = 0.7, sigmaT = 0.7, pve_x = 0, pve_g = 0.0001,
                                  maf, kinsetting, Genotype_each=NULL, pre_vars_step2, GRM) {
  
  list2env(pre_vars_step2, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype and covariate.
  if(!is.na(maf)){
    G_all <- generate_gmatrix_kinship(n_obs,miss,maf,kinsetting)
  }
  if(!is.null(Genotype_each)){
    G_all <- Genotype_each
  }
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all,X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - pve_g - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- sqrt(pve_g * V_total /var(G_all))
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  # 生成hatY
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  # 生成Y
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)
  
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
  # n_families: 家庭数量(5000个家庭=10000个人)
  # n_snps: SNP数量
  # maf: 次要等位基因频率
  # rho: 家庭成员间的基因型相关性
  
  # 1. 生成独立的标准正态分布变量(家庭间独立)
  # 每个SNP生成2n_families个独立样本
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  n_families <- n_all/2
  Z <- matrix(rnorm(n_snps * 2 * n_families), nrow = n_snps)
  
  # 2. 创建家庭成员间的相关性
  # 对于每个家庭中的两个人，他们的基因型是Z1和rho*Z1 + sqrt(1-rho^2)*Z2
  # 其中Z1和Z2是独立的
  Z_family <- matrix(0, nrow = n_snps, ncol = 2 * n_families)
  
  for (i in 1:n_families) {
    # 家庭中第一个人的基因型(完全独立)
    Z_family[, 2*i - 1] <- Z[, 2*i - 1]
    
    # 家庭中第二个人的基因型(与第一个人相关)
    Z_family[, 2*i] <- rho * Z[, 2*i - 1] + sqrt(1 - rho^2) * Z[, 2*i]
  }
  
  # 3. 转换为二项分布以控制MAF
  # 计算分位数阈值
  threshold <- qnorm(1 - maf)
  
  # 初始化基因型矩阵(0,1,2)
  genotypes <- matrix(0, nrow = 2 * n_families, ncol = n_snps)
  
  # 转换为基因型
  for (i in 1:n_snps) {
    # 对于每个SNP，将正态分布转换为基因型
    # 0: 低于阈值(主要等位基因纯合子)
    # 1: 在阈值和对称阈值之间(杂合子)
    # 2: 高于对称阈值(次要等位基因纯合子)
    
    # 使用对称的双侧阈值
    upper_threshold <- threshold
    lower_threshold <- -threshold
    
    # 转换为基因型
    gt <- ifelse(Z_family[i,] < lower_threshold, 0, 
                 ifelse(Z_family[i,] > upper_threshold, 2, 1))
    
    genotypes[, i] <- gt
  }
  
  # 4. 标准化基因型(均值为0，方差为1)
  # 首先中心化
  genotypes_centered <- scale(genotypes, center = TRUE, scale = FALSE)
  
  # 然后标准化方差
  genotypes_scaled <- genotypes_centered / sqrt(2 * maf * (1 - maf))
  
  # 转置矩阵，使行为个体，列为SNP
  genotypes_final <- t(genotypes_scaled)
  
  return(t(genotypes_final))
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
DGP_GRMr_gmodel_step3_rho <- function(n_obs,miss,tauT = 0.7, sigmaT = 0.7, pve_x = 0, pve_g = 0.0001,
                                      Genotype_each, pre_vars_step2, GRM) {
  
  list2env(pre_vars_step2, envir = environment())
  
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  
  # Genotype and covariate.
  G_all <- Genotype_each
  X_all <- stats::rnorm(n = n_all)
  Z_all <- cbind(G_all,X_all)
  
  # 计算随机效应和误差项的总变异
  r <- tauT^2 + sigmaT^2
  # 计算总表型变异
  V_total <- r / (1 - pve_g - pve_x)
  
  # 根据期望的贡献比例设置固定效应系数
  beta_G <- sqrt(pve_g * V_total /var(G_all))
  beta_X <- sqrt(pve_x * V_total/var(X_all)) 
  
  # 固定效应部分
  mu_all <- Z_all %*% c(beta_G, beta_X)
  
  # 生成hatY
  S <- as.numeric(mu_all + t(chol_Sigma22) %*% rnorm(n_all))
  
  # 生成Y
  mu_cond <- mu_all[1:n_obs] + Sigma12_Sigma22inv %*% (S - mu_all)
  epsl <- rnorm(n_all)
  Y_obs <- rep(NA,n_all)
  Y_obs[1:n_obs] <- as.numeric(mu_cond + t(L_cond) %*% epsl[1:n_obs])
  
  ### generate oracle Y 
  mu_cond_oracle <- mu_all + Sigma12_oracle_Sigma22inv %*% (S - mu_all)
  Y <- as.numeric(mu_cond_oracle + t(L_cond_oracle) %*% epsl)
  
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
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )%>% na.omit()
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df$X
  )
  
  model_lm1 <- lm(y ~ x, data = reml_data_SynSurrG)
  model_lm2 <- lm(haty ~ x, data = lm_data_SynSurrG, na.action = na.exclude)
  
  ## residuals aligned to their index sets. na.action = na.exclude keeps the NA
  ## padding, so res_S[obs_protein_index] refers to the right individuals even
  ## if lm() dropped rows; the old code indexed a shortened residual vector by
  ## position, which silently misaligns whenever a covariate has an NA.
  res_T <- as.numeric(model_lm1$residuals)
  res_S <- as.numeric(residuals(model_lm2))
  if (length(res_T) != n_obs)
    stop("score_test_SynSurrG_power_sim_one_g_to_one_Y: lm1 residual length ", length(res_T), " != n_obs ", n_obs,
         " -- NA in covariates or surrogate at observed individuals")
  if (length(res_S) != n)
    stop("score_test_SynSurrG_power_sim_one_g_to_one_Y: lm2 residual length ", length(res_S), " != n ", n)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(res_T) %*% GRM_o %*% res_T)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "score_test_SynSurrG_power_sim_one_g_to_one_Y")
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,
    dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  ## inv_Sigma11 is never used downstream (it is not in `params`): dropping it
  ## saves one full block-wise inversion per protein and removes a crash site.
  ## inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- matrix_inv_Amatrix(Amatrix = A)
  
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  #b3 <- Sigma12_Sigma22inv %*% hatY
  
  params <- c(
    "obs_protein_index","V11", 
    "Y", "hatY", "X_obs", "Att", "A22", "Atb", 
    "B1", "B_mat", "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11", "Sigma12_Sigma22inv", "Sigma12",
    "vc"
  )
  
  SynSurrG_step1_pars <- mget(params, envir = environment())
  
  G_all <- mydf$G_all
  
  results <- score_test_SynSurrG_single(G_all,step1_pars = SynSurrG_step1_pars)
  
  attr(results, "vc") <- vc
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
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  n_obs <- length(obs_protein_index)
  
  data_obsG <- data.frame(
    y = test_df$Y_obs,
    x = test_df$X
  )%>% na.omit()
  
  model_lm1 <- lm(y ~ x, data = data_obsG)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm1$residuals) %*% GRM_o %*% model_lm1$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm1$residuals^2), label = "score_test_ObsG_power_sim_one_g_to_one_Y")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # 计算 Σ11^{-1}
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  
  X_all <- data.frame(intercept = rep(1,n),test_df$X)
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  G_all <- mydf$G_all
  
  SX <- inv_Sigma11 %*% X_obs
  XtSX_inv <- solve(crossprod(X_obs, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X_obs %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  ## for obs test statistic
  G_obs <- G_all[obs_protein_index]
  SU <- as.numeric(crossprod(G_obs, invS11_res))
  A <- crossprod(SX, G_obs)
  VU <- as.numeric(colSums(G_obs * (inv_Sigma11 %*% G_obs)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_ObsG","negative_log10_pval_ObsG")
  
  attr(results, "vc") <- vc
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
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  GRM <- mydf$GRM
  
  lm_data_OracleG <- data.frame(
    y = test_df$Y_all,
    x = test_df$X
  )
  
  model_lm <- lm(y ~ x, data = lm_data_OracleG)
  
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_oall %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "score_test_OracleG_power_sim_one_g_to_one_Y")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM + hat_sigma_T2 * Diagonal(n)
  
  inv_Sigma11 <- matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- lm_data_OracleG$y
  X <- cbind(rep(1,n),lm_data_OracleG$x)
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  beta_hat <- XtSX_inv %*% (t(SX) %*% Y)
  residual <- Y - X %*% beta_hat
  invS11_res <- inv_Sigma11 %*% residual
  
  
  ## G_all was never assigned here (the ObsG twin does `G_all <- mydf$G_all`),
  ## so this function silently depended on a G_all leaking in from the caller.
  G_all <- mydf$G_all
  
  SU <- as.numeric(crossprod(G_all, invS11_res))
  A <- crossprod(SX, G_all)
  VU <- as.numeric(colSums(G_all * (inv_Sigma11 %*% G_all)) - colSums(A * (XtSX_inv %*% A)))
  
  results <- compute_score(SU, VU)
  
  names(results) <- c("T_score_OracleG","negative_log10_pval_OracleG")
  
  attr(results, "vc") <- vc
  return(results)
}


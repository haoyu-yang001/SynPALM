## functions ########
library(lattice)
library(data.table)
library(dplyr)
library(Matrix)
library(GMMAT)
library(ranger)
library(parallel)
library(withr)

## ###########################################################################
## funcs_v3.R -- funcs_v2.R plus a numerical-stability layer (2026-08-22)
##
## Verified: on data where no constraint binds, every function below
## reproduces the funcs_v2.R output to ~3e-15, so existing results do not
## move. What changes is only the cases where v2 crashed or silently produced
## an indefinite matrix.
##
##   options(synsurrg.on_nonpsd = "project")  # "project" (default) | "stop"
##   options(synsurrg.verbose   = TRUE)
##   options(synsurrg.ridge_e   = 0.01)       # residual-variance floor, as a
##                                            # fraction of mean(resid^2)
##   options(synsurrg.ncores    = NULL)       # NULL -> SLURM_CPUS_PER_TASK, else 2
##
## See the block comment above .ss_opt() for the full rationale.
## ###########################################################################


reset_temp_dir <- function(dir) {
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE) # 删除现有目录
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE) # 重新创建
}

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
  
  message("✅ .sample file written to: ", file_out)
}

# hh: list，每个元素是一个 block 的样本下标向量（如你给的 hh[[1]], hh[[2]] ...）
# prop: 训练集比例（默认 0.8）
# iterations: 随机搜索次数，越大越接近目标（1e3~1e4 一般就够）
# seed: 可选，设定随机种子以复现
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

merge_results <- function(res_list) {
  nm <- names(res_list[[1]])
  setNames(
    lapply(nm, function(x) do.call(c, lapply(res_list, `[[`, x))),
    nm
  )
}

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


## ===========================================================================
##  Numerical-stability layer  (added 2026-08-22)
##
##  Two problems this fixes.
##
##  (1) matrix_inv_block() / matrix_inv_Amatrix() had no error handling. When a
##      block failed to invert, mclapply() returned a "try-error" object in that
##      list position and bdiag() then died with
##          no method or default for coercing "try-error" to "CsparseMatrix"
##      which hides the real cause. Errors are now caught per block and
##      reported with the block id, its size and its smallest eigenvalue;
##      non-PD blocks get an escalating ridge, then an eigenvalue-floored
##      pseudo-inverse.
##
##  (2) The Haseman-Elston moment estimators of the variance components are
##      unconstrained, so tau^2 and sigma^2 can come out negative. Per
##      individual the joint covariance of (T_i, S_i) is
##
##            G (x) GRM  +  E (x) I ,
##
##            G = [ tau_T2   tau_TS  ]      E = [ sigma_T2   sigma_TS  ]
##                [ tau_TS   tau_S2  ]          [ sigma_TS   sigma_S2  ]
##
##      and because GRM >= 0 this is positive definite iff G >= 0 and E > 0,
##      i.e. tau_TS^2 <= tau_T2 * tau_S2 (and the same for the sigmas).
##      Flooring the six numbers one at a time with max(., 0.001) does NOT
##      give that: tau_T2 = 0.30, tau_S2 = 0.001, tau_TS = 0.25 passes every
##      floor yet has det(G) < 0. And since
##
##            A = Sigma11 - Sigma12 Sigma22^{-1} Sigma21
##
##      is the Schur complement of the joint matrix, A > 0 iff the joint is
##      > 0 -- so an indefinite E does not show up in Sigma11 or Sigma22, it
##      shows up later in matrix_inv_Amatrix(). The 2x2 matrices are now
##      projected onto the PSD cone instead, which guarantees Sigma11,
##      Sigma22 AND A are all positive definite.
##
##  Behaviour is controlled by options(), so no call site had to change:
##
##      options(synsurrg.on_nonpsd = "project")  # "project" (default) | "stop" | "none"
##                                               # "none" = pre-2026-08-22 behaviour:
##                                               # raw estimates, only sigma_T2
##                                               # floored at 0.01 (added 2026-10-08)
##      options(synsurrg.verbose   = TRUE)       # print VC / block diagnostics
##      options(synsurrg.ridge_e   = 0.01)       # residual-variance floor, as a
##                                               # fraction of mean(resid^2)
##                                               # (replaces the old absolute
##                                               #  max(., 0.01) / max(., 0.001))
##      options(synsurrg.ncores    = NULL)       # NULL -> SLURM_CPUS_PER_TASK, else 2
##
##  Every *_step1 / *_estimate function now also returns a `vc` element holding
##  the raw and the used variance components plus the projection flags, so a
##  regularised protein is identifiable afterwards instead of silent.
## ===========================================================================

.ss_opt <- function(name, default) {
  v <- getOption(paste0("synsurrg.", name), default)
  if (is.null(v)) default else v
}

.ss_cores <- function() {
  v <- getOption("synsurrg.ncores", NULL)
  if (!is.null(v) && is.finite(v) && v >= 1) return(as.integer(v))
  x <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (!is.na(x) && x >= 1L) return(x)
  2L
}

.ss_verbose <- function() isTRUE(.ss_opt("verbose", TRUE))
.ss_strict  <- function() identical(.ss_opt("on_nonpsd", "project"), "stop")

## ---- method-of-moments ratio ------------------------------------------------
##  tau = a1 / sum(K_offdiag^2). Without any related pair among the individuals
##  used (e.g. no relatives among the labelled) the denominator is 0 and the
##  genetic component cannot be estimated from relatives: it is set to 0, which
##  reduces the model to the bivariate linear model of SynSurr (2026-10-13).
.ss_mom_ratio <- function(a1, denom) {
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    if (.ss_verbose()) message("no related pairs among the individuals used: genetic variance component set to 0")
    return(0)
  }
  a1 / denom
}

## ---- project a symmetric 2x2 onto the PSD cone -----------------------------
##  floor_diag: lower bound applied to the diagonal AFTER the projection
##              (adding a non-negative diagonal preserves PSD).
##  shrink    : pulls the correlation just inside the unit circle so the result
##              is strictly PD rather than only PSD.
##  Both are no-ops for a matrix that is already comfortably positive definite,
##  so a well-behaved protein gets numerically identical estimates to before.
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

## ---- bivariate (T, S) variance components ----------------------------------
##  var_T / var_S: mean(resid^2) for the two traits, used to scale the ridge so
##  the floor is relative to the trait scale rather than a hard-coded 0.01.
constrain_vc_bivariate <- function(tau_T2, tau_TS, tau_S2,
                                   sigma_T2, sigma_TS, sigma_S2,
                                   var_T = 1, var_S = 1, label = "",
                                   dbar_T = 1, dbar_S = 1) {
  raw <- c(tau_T2 = tau_T2, tau_TS = tau_TS, tau_S2 = tau_S2,
           sigma_T2 = sigma_T2, sigma_TS = sigma_TS, sigma_S2 = sigma_S2)
  
  eps   <- .ss_opt("ridge_e", 0.01)
  scale <- mean(c(abs(var_T), abs(var_S)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  
  g <- psd_project_2x2(tau_T2,   tau_TS,   tau_S2,   floor_diag = 0)
  ## Re-solve the residual components from the same diagonal moment equations
  ## with the constrained genetic block (2026-10-13):
  ##   a2/n = tau * dbar + sigma  =>  sigma = sigma_raw + (tau_raw - tau_used) * dbar
  ## with dbar = tr(GRM)/n of the matrix used in that equation (T and TS: the
  ## observed block, S: all individuals). Keeping sigma_raw instead inflates the
  ## residual variance by |tau_raw| * dbar whenever tau_raw < 0, which makes the
  ## score test conservative (lambda_GC 0.85-0.93 under permutation).
  ## A no-op when g is unchanged. options(synsurrg.rematch_sigma = FALSE)
  ## restores the 2026-08-22 behaviour.
  if (isTRUE(.ss_opt("rematch_sigma", TRUE)) && g$changed) {
    sigma_T2 <- sigma_T2 + (tau_T2 - g$a)  * dbar_T
    sigma_TS <- sigma_TS + (tau_TS - g$ab) * dbar_T
    sigma_S2 <- sigma_S2 + (tau_S2 - g$b)  * dbar_S
  }
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
             rematch_sigma = isTRUE(.ss_opt("rematch_sigma", TRUE)),
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

## ---- univariate (T only) variance components -------------------------------
##  Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0,
##  so here a simple floor is the correct constraint (no Cauchy-Schwarz to
##  worry about).
constrain_vc_univariate <- function(tau2, sigma2, var_T = 1, label = "", dbar = 1) {
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
  
  ## tau2 floored at 0; when it is, sigma2 is re-solved from
  ## a2/n = tau2 * dbar + sigma2 (see constrain_vc_bivariate; 2026-10-13).
  ## sigma2 is then floored as before (tau2 is left alone, as in the bivariate case).
  rematch   <- isTRUE(.ss_opt("rematch_sigma", TRUE)) && is.finite(tau2) && is.finite(sigma2)
  bad_tau   <- !is.finite(tau2) || tau2 < 0
  tau_used  <- if (bad_tau) 0 else tau2
  sig_used  <- if (rematch) sigma2 + (tau2 - tau_used) * dbar else sigma2
  bad_sigma <- !is.finite(sig_used) || sig_used < eps * scale
  if (bad_sigma) sig_used <- eps * scale
  tau2 <- tau_used; sigma2 <- sig_used
  
  if ((bad_tau || bad_sigma) && .ss_strict())
    stop(sprintf("[%s] variance components out of range: %s", label,
                 paste(sprintf("%s=%+.4g", names(raw), raw), collapse = " ")),
         call. = FALSE)
  
  if (.ss_verbose() && (bad_tau || bad_sigma))
    message(sprintf("[%s] VC floored: tau_T2 %+.4g -> %.4g, sigma_T2 %+.4g -> %.4g",
                    label, raw[["tau_T2"]], tau2, raw[["sigma_T2"]], sigma2))
  
  list(tau2 = tau2, sigma2 = sigma2,
       vc = list(raw = raw, used = c(tau_T2 = tau2, sigma_T2 = sigma2),
                 floored_tau = bad_tau, floored_sigma = bad_sigma,
                 rematch_sigma = rematch, label = label))
}

## ---- greedy merge of consecutive blocks ------------------------------------
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

## ---- invert ONE block: chol -> chol+ridge -> eigenvalue-floored inverse ----
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

## ---- parallel driver: invert all blocks, report, assemble ------------------
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

find_blocks_vectorized_threshold <- function(A, threshold = 0.005) {
  find_blocks_vectorized(A, threshold = threshold)
}

block_matrix <- function(tl, tr, bl, br) {
  rbind(cbind(tl, tr), cbind(bl, br))
}

matrix_inv_block <- function(wait_matrix, threshold = 1000, label = "Sigma") {
  blocks        <- find_blocks_vectorized(wait_matrix)
  merged_blocks <- .ss_merge_blocks(blocks, threshold)
  .ss_block_inverse(wait_matrix, merged_blocks, label = label)
}

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

compute_score <- function(Ug, Vu) {
  score <- as.numeric(Ug^2 / Vu)
  negative_log10_pval <- -pchisq(score, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  list(score = score, negative_log10_pval = negative_log10_pval)
}

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

#################### functions for simulation #################
## Genetic and residual covariances of (Y, S) for a target correlation rho,
## split 3/5 genetic : 2/5 residual as in the submitted simulations. When the
## genetic share would break Cauchy-Schwarz (|tauTS| > tauT * tauS; rho = 0.75
## with the default components gives 0.285 > 0.280), the genetic covariance is
## set to `cap` x the bound and the remainder goes to the residual covariance,
## so corr(Y, S) is still rho (Reviewer 2, comment 5; 2026-10-13).
## tauT, tauS, sigmaT, sigmaS are standard deviations, as everywhere below.
sim_cov_TS <- function(rho, tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5,
                       g_share = 3/5, cap = 0.99) {
  total   <- sqrt((tauT^2 + sigmaT^2) * (tauS^2 + sigmaS^2)) * rho
  tauTS   <- total * g_share
  bound   <- cap * tauT * tauS
  if (abs(tauTS) > bound) tauTS <- sign(tauTS) * bound
  sigmaTS <- total - tauTS
  if (abs(sigmaTS) >= sigmaT * sigmaS)
    stop("sim_cov_TS: rho = ", rho, " cannot be reached with these variance components")
  c(tauTS = tauTS, sigmaTS = sigmaTS,
    rho_g = tauTS / (tauT * tauS), rho_e = sigmaTS / (sigmaT * sigmaS))
}
#### for type 1 error control #########

##目前的这种grm是两两0.5的结构
DGP <- function(n_obs,miss,rho,maf = 0.25,tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5, 
                pve_g = 0, pve_x = 0.10) {
  
  ## using rho to generate the tauTS and sigmaTS
  cv <- sim_cov_TS(rho, tauT, tauS, sigmaT, sigmaS)
  tauTS <- cv[["tauTS"]]; sigmaTS <- cv[["sigmaTS"]]
  
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
  cv <- sim_cov_TS(rho, tauT, tauS, sigmaT, sigmaS)
  tauTS <- cv[["tauTS"]]; sigmaTS <- cv[["sigmaTS"]]
  
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

typeIerror_generate_gmatrix <- function(n_obs,miss,chunk_size,maf){
  n_miss <- round(n_obs * miss / (1 - miss))
  n_all <- n_obs + n_miss
  g_matrix <- matrix(rbinom(n_all * chunk_size, size = 2, prob = maf), nrow = n_all, ncol = chunk_size)
  return(g_matrix)
}

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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_typeIerror_step1",
                                dbar_T = mean(diag(GRM_obs_obs)), dbar_S = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_ablation_estimate",
                                dbar_T = mean(diag(GRM_obs_obs)), dbar_S = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "OracleG_typeIerror_step1",
                                                                     dbar = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "OracleG_ablation_estimate",
                                                                     dbar = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_typeIerror_step1",
                                                                     dbar = mean(diag(GRM_obs_obs)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_hatY_typeIerror_step1",
                                                                     dbar = mean(diag(GRM_obs_obs)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_uno)^2)) #equal to sum(diag((GRM_uno) %*% t(GRM_unobs_unobs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_unobs_unobs)))/nrow(GRM_unobs_unobs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "unObsG_hatY_typeIerror_step1",
                                                                     dbar = mean(diag(GRM_unobs_unobs)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_eachG_step1",
                                dbar_T = mean(diag(GRM_obs_obs)), dbar_S = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_ablation_estimate",
                                                                     dbar = mean(diag(GRM_obs_obs)))
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

############# use independent samples ###########

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

#### for power compare #########

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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "score_test_SynSurrG_power_sim_one_g_to_one_Y",
                                dbar_T = mean(diag(GRM_obs_obs)), dbar_S = mean(diag(GRM)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_o)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm1$residuals * model_lm1$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm1$residuals^2), label = "score_test_ObsG_power_sim_one_g_to_one_Y",
                                                                      dbar = mean(diag(GRM_obs_obs)))
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
  hat_tau_T2 <- .ss_mom_ratio(a1, sum((GRM_oall)^2)) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM)))/nrow(GRM)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "score_test_OracleG_power_sim_one_g_to_one_Y",
                                                                     dbar = mean(diag(GRM)))
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
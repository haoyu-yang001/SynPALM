## ---------------------------------------------------------------------------
## SynPALM end-to-end pipeline.
##
##   synpalm_folds()    relatedness-cluster-based K-fold assignment
##   synpalm_predict()  cross-fitted random-forest synthetic phenotype
##   synpalm_null()     null model / variance components (LMM step 1)
##   synpalm_scan()     genome-wide score tests (LMM step 2)
##   synpalm_gwas()     all of the above in one call
##
## Each step is a port of the analysis scripts used for the UK Biobank results
## (make_folds_once.R, step1_protein_crossfit.R, step2_scan.R); the model
## fitting and the score tests are the functions in synpalm_functions.R.
## ---------------------------------------------------------------------------

## Evaluate `expr` with a fixed seed and restore the caller's RNG state.
.with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

## Sample IDs of an object: names() for vectors, rownames() otherwise.
.ids_of <- function(x, what) {
  id <- if (is.null(dim(x))) names(x) else rownames(x)
  if (is.null(id)) stop(what, " must carry sample IDs as ",
                        if (is.null(dim(x))) "names" else "rownames", call. = FALSE)
  as.character(id)
}

.check_grm <- function(grm) {
  if (!inherits(grm, "Matrix")) grm <- Matrix::Matrix(grm, sparse = TRUE)
  if (nrow(grm) != ncol(grm)) stop("grm must be square", call. = FALSE)
  if (is.null(rownames(grm))) stop("grm must carry sample IDs as rownames", call. = FALSE)
  if (is(grm, "symmetricMatrix")) grm <- as(grm, "generalMatrix")
  grm
}

## Connected components of the GRM graph, labelled by the smallest index in
## each component. Label propagation with pointer jumping; relatedness
## clusters have a small diameter, so this converges in a few sweeps.
.grm_components <- function(grm) {
  sm <- summary(grm)
  i <- sm$i; j <- sm$j
  lab <- seq_len(nrow(grm))
  repeat {
    m <- pmin(lab[i], lab[j])
    o <- order(m, decreasing = TRUE)        # last write wins: the smallest label
    new <- lab
    new[i[o]] <- pmin(lab[i[o]], m[o])
    new <- new[new]
    if (identical(new, lab)) break
    lab <- new
  }
  lab
}

## ===========================================================================
## Folds
## ===========================================================================

#' Relatedness-Aware Cross-Fitting Folds
#'
#' Splits the analysis cohort into \code{K} folds so that every relatedness
#' cluster (connected component of the sparse GRM) lies entirely within one
#' fold. Multi-member clusters are placed largest-first into the emptiest
#' fold, then singletons are dealt out to equalise fold sizes. Because whole
#' clusters are assigned, the out-of-fold prediction for an individual never
#' uses the phenotype of any of their relatives.
#'
#' One individual per cluster is also drawn as the set of independent
#' individuals used by the unrelated-sample methods and by
#' \code{\link{fix_constant_columns}}.
#'
#' The block-wise algorithms require every cluster to occupy a contiguous
#' range of rows. The individuals are therefore reordered so that each
#' cluster is contiguous; an order that already has this property (such as
#' the UK Biobank sparse GRM) is left unchanged. The returned \code{ids} are
#' the analysis order used by every later step.
#'
#' @param grm Sparse genetic relatedness matrix with sample IDs as dimnames.
#' @param K Number of folds.
#' @param seed Seed for the singleton assignment.
#' @param independent_seed Seed for drawing one individual per cluster.
#'
#' @return A list of class \code{"synpalm_folds"} with \code{ids} (the
#'   analysis order), \code{fold} (integer, aligned to \code{ids}),
#'   \code{blocks} (relatedness clusters as index vectors into \code{ids}),
#'   \code{independent_indices}, \code{K} and \code{seed}.
#' @export
synpalm_folds <- function(grm, K = 5L, seed = 20260813L, independent_seed = 25324L) {
  grm  <- .check_grm(grm)
  comp <- .grm_components(grm)
  ord  <- order(comp, seq_along(comp))
  if (is.unsorted(ord)) grm <- grm[ord, ord]
  blocks <- find_blocks_vectorized(grm)
  if (length(blocks) != length(unique(comp)))
    stop("internal error: block detection disagrees with the GRM components", call. = FALSE)
  sz     <- lengths(blocks)
  n_all  <- nrow(grm)

  fold_of_block <- integer(length(blocks))
  load_k        <- numeric(K)

  multi <- which(sz > 1)
  for (b in multi[order(sz[multi], decreasing = TRUE)]) {
    k <- which.min(load_k)
    fold_of_block[b] <- k
    load_k[k] <- load_k[k] + sz[b]
  }

  sing <- which(sz == 1)
  n_s  <- length(sing)
  .with_seed(seed, {
    perm <- sample(sing)
  })
  target <- sum(sz) / K
  need   <- as.integer(floor(pmax(0, target - load_k)))
  if (sum(need) > n_s) need <- as.integer(round(need * n_s / sum(need)))
  need[K] <- n_s - sum(need[-K])
  if (any(need < 0)) {
    fold_of_block[perm] <- rep_len(seq_len(K), n_s)
  } else {
    fold_of_block[perm] <- rep(seq_len(K), times = need)
  }

  fold_vec <- integer(n_all)
  for (b in seq_along(blocks)) fold_vec[blocks[[b]]] <- fold_of_block[b]
  stopifnot(all(fold_vec %in% seq_len(K)))
  split_clusters <- sum(vapply(blocks, function(ix) length(unique(fold_vec[ix])), integer(1)) > 1)
  if (split_clusters > 0) stop("internal error: ", split_clusters, " clusters split across folds")

  independent_indices <- .with_seed(independent_seed, {
    vapply(blocks, function(x) x[sample.int(length(x), 1)], integer(1))
  })

  structure(list(ids = rownames(grm), fold = fold_vec, blocks = blocks,
                 independent_indices = independent_indices, K = K, seed = seed),
            class = "synpalm_folds")
}

## ===========================================================================
## Step 1: cross-fitted synthetic phenotype
## ===========================================================================

#' Cross-Fitted Random-Forest Synthetic Phenotype
#'
#' Predicts the protein for every individual with a random forest trained,
#' within each fold, only on labelled individuals outside that fold. The whole
#' supervised pipeline, feature screening included, is refitted inside each
#' fold, so no individual's prediction depends on their own phenotype or on
#' that of any relative.
#'
#' Screening, per fold and on the training individuals only:
#' \itemize{
#'   \item columns of \code{rf_features} other than \code{rf_fixed}: the
#'     \code{n_top} with the largest absolute Pearson correlation with the
#'     protein; the \code{rf_fixed} columns are then always added;
#'   \item columns of \code{rf_binary} (0/1 indicators such as disease
#'     categories): ranked by two-sided Wilcoxon \eqn{-\log_{10} p}; the top
#'     \code{binary_top[1]} are kept if at most \code{binary_top[1]} exceed
#'     \code{binary_sig}, otherwise the top \code{binary_top[2]}. They enter
#'     the forest as factors.
#' }
#' With \code{screen = FALSE} every column is used.
#'
#' @param protein Numeric vector of the measured protein, \code{NA} where
#'   unmeasured, named by sample ID.
#' @param rf_features Data frame of candidate predictors, rownames = sample IDs.
#' @param folds Output of \code{\link{synpalm_folds}}. Its \code{ids} define
#'   the analysis set and order.
#' @param rf_binary Optional data frame of 0/1 predictors, rownames = sample IDs.
#' @param rf_fixed Names of \code{rf_features} columns that are always kept
#'   and are not screening candidates (e.g. age and sex).
#' @param screen Logical; screen predictors within each fold.
#' @param n_top Number of \code{rf_features} columns kept by screening.
#' @param binary_sig,binary_top Screening rule for \code{rf_binary}, see Details.
#' @param num.trees Number of trees per forest.
#' @param seed Fold \code{k} uses \code{seed + k}.
#' @param num.threads Threads for \code{ranger}.
#' @param verbose Print per-fold progress.
#'
#' @return A list of class \code{"synpalm_prediction"} with \code{ids},
#'   \code{protein} (aligned to \code{ids}), \code{protein_hat},
#'   \code{rho_oof} (Pearson correlation between the synthetic and the measured
#'   protein over all labelled individuals, each predicted out of fold),
#'   \code{rho_by_fold}, \code{accuracy}, \code{n_obs} and \code{fold_features}.
#'   \code{accuracy} is a data frame with one row per fold plus a row
#'   \code{fold = "all"}: \code{n_train} (labelled individuals the fold's forest
#'   was trained on), \code{n_labelled} (labelled individuals in the fold),
#'   \code{n_predicted} (all individuals in the fold), \code{rho},
#'   \code{rho_spearman} and \code{r2 = rho^2}.
#' @export
synpalm_predict <- function(protein, rf_features, folds, rf_binary = NULL, rf_fixed = NULL,
                            screen = TRUE, n_top = 100L, binary_sig = 5, binary_top = c(5L, 10L),
                            num.trees = 300L, seed = 260611L, num.threads = 1L, verbose = TRUE) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("synpalm_predict() needs the 'ranger' package", call. = FALSE)
  stopifnot(inherits(folds, "synpalm_folds"))
  ids <- folds$ids
  K   <- folds$K
  fold_vec <- folds$fold
  n_all    <- length(ids)

  take <- function(x, what) {
    m <- match(ids, .ids_of(x, what))
    if (anyNA(m)) stop(sum(is.na(m)), " analysis individuals are missing from ", what, call. = FALSE)
    if (is.null(dim(x))) x[m] else x[m, , drop = FALSE]
  }
  y_all <- as.numeric(take(protein, "protein"))
  feat  <- as.data.frame(take(rf_features, "rf_features"))
  dis   <- if (is.null(rf_binary)) NULL else as.data.frame(take(rf_binary, "rf_binary"))
  if (!all(rf_fixed %in% names(feat)))
    stop("rf_fixed not found in rf_features: ",
         paste(setdiff(rf_fixed, names(feat)), collapse = ", "), call. = FALSE)

  obs_index <- which(!is.na(y_all))
  if (length(obs_index) < 10) stop("fewer than 10 individuals with a measured protein", call. = FALSE)

  X_obs <- as.matrix(feat[obs_index, setdiff(names(feat), rf_fixed), drop = FALSE])
  Z_obs <- if (is.null(dis)) NULL else as.matrix(dis[obs_index, , drop = FALSE])
  y_obs <- y_all[obs_index]

  sanitize <- function(nm) gsub("[- ]", "_", nm)

  screen_fold <- function(train_rows) {
    if (!screen)
      return(list(cov_pred = names(feat), cov_dis = if (is.null(dis)) character(0) else names(dis),
                  n_sig = NA_integer_))
    tp <- match(train_rows, obs_index)

    r <- suppressWarnings(stats::cor(X_obs[tp, , drop = FALSE], y_obs[tp],
                                     use = "pairwise.complete.obs"))
    r[is.na(r)] <- 0
    cov_pred <- c(colnames(X_obs)[order(abs(r), decreasing = TRUE)[seq_len(min(n_top, ncol(X_obs)))]],
                  rf_fixed)

    if (is.null(Z_obs)) return(list(cov_pred = cov_pred, cov_dis = character(0), n_sig = NA_integer_))
    yv  <- y_obs[tp]
    sig <- vapply(seq_len(ncol(Z_obs)), function(i) {
      zi   <- Z_obs[tp, i]
      vals <- unique(zi[!is.na(zi)])
      if (length(vals) == 2) {
        p <- stats::wilcox.test(x = yv[zi == vals[2]], y = yv[zi == vals[1]],
                                alternative = "two.sided", paired = FALSE, exact = FALSE)$p.value
        -log10(p)
      } else 0
    }, numeric(1))
    n_sig   <- sum(sig > binary_sig)
    n_keep  <- if (n_sig <= binary_top[1]) binary_top[1] else binary_top[2]
    cov_dis <- colnames(Z_obs)[order(sig, decreasing = TRUE)[seq_len(min(n_keep, ncol(Z_obs)))]]
    list(cov_pred = cov_pred, cov_dis = cov_dis, n_sig = n_sig)
  }

  build_frame <- function(rows, cov_pred, cov_dis) {
    out <- feat[rows, cov_pred, drop = FALSE]
    if (length(cov_dis)) {
      b <- dis[rows, cov_dis, drop = FALSE]
      b[] <- lapply(b, function(x) factor(x, levels = c(0, 1), labels = c("0", "1")))
      out <- cbind(out, b)
    }
    names(out) <- sanitize(names(out))
    out
  }

  protein_hat   <- rep(NA_real_, n_all)
  fold_features <- vector("list", K)

  for (k in seq_len(K)) {
    train_rows <- intersect(obs_index, which(fold_vec != k))
    pred_rows  <- which(fold_vec == k)
    if (verbose) message(sprintf("fold %d: train on %d labelled, predict %d individuals",
                                 k, length(train_rows), length(pred_rows)))

    sel <- screen_fold(train_rows)
    fold_features[[k]] <- sel

    tr <- build_frame(train_rows, sel$cov_pred, sel$cov_dis)
    tr$protein <- y_all[train_rows]

    keep <- vapply(tr[setdiff(names(tr), "protein")],
                   function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
    drop_nm <- names(keep)[!keep]
    if (length(drop_nm)) {
      if (verbose) message("  dropping non-varying predictors: ", paste(drop_nm, collapse = ", "))
      tr <- tr[, c(setdiff(names(tr), c(drop_nm, "protein")), "protein"), drop = FALSE]
    }

    rf <- ranger::ranger(protein ~ ., data = tr, num.trees = num.trees, importance = "impurity",
                         write.forest = TRUE, seed = seed + k, num.threads = num.threads)

    pr <- build_frame(pred_rows, sel$cov_pred, sel$cov_dis)
    pr <- pr[, setdiff(names(tr), "protein"), drop = FALSE]
    protein_hat[pred_rows] <- stats::predict(rf, data = pr, num.threads = num.threads)$predictions
    rm(tr, pr, rf)
  }
  stopifnot(!anyNA(protein_hat))

  ## accuracy of the synthetic phenotype: correlation with the measured protein
  ## among labelled individuals, each predicted by a forest that never saw them
  acc_row <- function(fold, ii, n_train, n_pred) {
    r <- if (length(ii) > 2) stats::cor(y_all[ii], protein_hat[ii]) else NA_real_
    data.frame(fold = fold, n_train = n_train, n_labelled = length(ii), n_predicted = n_pred,
               rho = r, rho_spearman = if (length(ii) > 2)
                 stats::cor(y_all[ii], protein_hat[ii], method = "spearman") else NA_real_,
               r2 = r^2, stringsAsFactors = FALSE)
  }
  accuracy <- do.call(rbind, c(
    lapply(seq_len(K), function(k)
      acc_row(as.character(k), intersect(obs_index, which(fold_vec == k)),
              length(intersect(obs_index, which(fold_vec != k))), sum(fold_vec == k))),
    list(acc_row("all", obs_index, NA_integer_, n_all))))
  rho_oof  <- accuracy$rho[K + 1]
  rho_fold <- accuracy$rho[seq_len(K)]
  if (verbose) message(sprintf("synthetic vs measured protein, %d-fold cross-fit: rho = %.4f (R2 = %.4f)",
                               K, rho_oof, rho_oof^2))

  structure(list(ids = ids, protein = y_all, protein_hat = protein_hat, rho_oof = rho_oof,
                 rho_by_fold = rho_fold, accuracy = accuracy, n_obs = length(obs_index),
                 fold_features = fold_features),
            class = "synpalm_prediction")
}

## ===========================================================================
## Step 2a: null model
## ===========================================================================

.vc_modes <- c("none", "project", "stop")

#' Fit the SynPALM Null Model
#'
#' Inverse-normal transforms the observed protein (observed individuals only)
#' and the synthetic phenotype (all individuals), estimates the variance
#' components and precomputes everything the score tests need.
#'
#' @param protein Numeric vector, \code{NA} where unmeasured, aligned to
#'   \code{ids} (e.g. \code{prediction$protein}).
#' @param synthetic Numeric vector of the synthetic phenotype, aligned to
#'   \code{ids} (e.g. \code{prediction$protein_hat}).
#' @param covariates Data frame of adjustment covariates (age, sex, PCs, ...),
#'   rownames = sample IDs. Factor and character columns are expanded into
#'   indicator columns. No missing values are allowed.
#' @param grm Sparse GRM with sample IDs as dimnames.
#' @param folds Output of \code{\link{synpalm_folds}}; supplies the analysis
#'   order and the independent individuals.
#' @param methods Any of \code{"SynSurrG"} (SynPALM), \code{"ObsG"} (LMM on
#'   observed individuals only), \code{"SynSurr"} and \code{"Obs"} (the same on
#'   one individual per relatedness cluster, without the GRM).
#' @param vc_constraint How invalid variance-component estimates are handled:
#'   \code{"none"} keeps the raw estimates (residual variance of the protein
#'   floored at 0.01), \code{"project"} projects them onto the PSD cone,
#'   \code{"stop"} raises an error. See \code{\link{SynPALM-package}}.
#'
#' @return A list of class \code{"synpalm_null"} with \code{ids},
#'   \code{obs_protein_index}, \code{independent_indices}, \code{methods},
#'   \code{pars} (one step-1 parameter list per method, each with a \code{vc}
#'   element where applicable) and \code{vc_constraint}.
#' @export
synpalm_null <- function(protein, synthetic, covariates, grm, folds,
                         methods = c("SynSurrG", "ObsG"), vc_constraint = "none") {
  stopifnot(inherits(folds, "synpalm_folds"))
  methods <- match.arg(methods, c("SynSurrG", "ObsG", "SynSurr", "Obs"), several.ok = TRUE)
  vc_constraint <- match.arg(vc_constraint, .vc_modes)
  ids <- folds$ids
  n   <- length(ids)
  if (length(protein) != n || length(synthetic) != n)
    stop("protein and synthetic must be aligned to folds$ids", call. = FALSE)

  grm <- .check_grm(grm)
  m <- match(ids, rownames(grm))
  if (anyNA(m)) stop(sum(is.na(m)), " analysis individuals are missing from grm", call. = FALSE)
  GRM <- grm[m, m]

  mc <- match(ids, .ids_of(covariates, "covariates"))
  if (anyNA(mc)) stop(sum(is.na(mc)), " analysis individuals are missing from covariates", call. = FALSE)
  cv <- as.data.frame(covariates)[mc, , drop = FALSE]
  if (anyNA(cv)) stop("covariates contain missing values for ",
                      sum(!stats::complete.cases(cv)), " analysis individuals", call. = FALSE)
  X <- stats::model.matrix(~ ., data = cv)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  colnames(X) <- make.names(gsub("-", "_", colnames(X)), unique = TRUE)

  ## transforms as in step2_scan.R
  ## INT() subsets rows with data[-i, ], which needs at least two columns
  tmp <- INT(data.frame(protein = as.numeric(protein), row = seq_len(n)), "protein")
  r   <- rank(as.numeric(synthetic))
  yhat_int <- stats::qnorm((r - 0.375) / (n - 2 * 0.375 + 1))

  mydf <- list(X_all = as.data.frame(X), S = yhat_int, Y_obs = tmp$int, GRM = GRM)
  independent_indices <- folds$independent_indices

  old <- options(synsurrg.on_nonpsd = vc_constraint)
  on.exit(options(old), add = TRUE)

  pars <- list()
  if ("SynSurrG" %in% methods) pars$SynSurrG <- SynSurrG_ablation_estimate(mydf)
  if ("ObsG"     %in% methods) pars$ObsG     <- ObsG_ablation_estimate(mydf)
  if ("SynSurr"  %in% methods) pars$SynSurr  <- SynSurr_ablation_estimate(mydf, independent_indices)
  if ("Obs"      %in% methods) pars$Obs      <- Obs_ablation_estimate(mydf, independent_indices)

  structure(list(ids = ids, obs_protein_index = which(!is.na(mydf$Y_obs)),
                 independent_indices = independent_indices, methods = methods,
                 pars = pars, vc_constraint = vc_constraint),
            class = "synpalm_null")
}

## ===========================================================================
## Step 2b: genome-wide scan
## ===========================================================================

.score_fun <- list(
  SynSurrG = function(G, p, ind) score_test_SynSurrG_multiply(G, p),
  ObsG     = function(G, p, ind) score_test_ObsG_multiply(G, p),
  SynSurr  = function(G, p, ind) score_test_SynSurr_multiply(G, p, ind),
  Obs      = function(G, p, ind) score_test_Obs_multiply(G, p, ind)
)

## Open a genotype source: a PLINK .bed path, or any matrix-like object with
## sample IDs as rownames and variant IDs as colnames.
.open_genotype <- function(genotype) {
  if (is.character(genotype) && length(genotype) == 1) {
    if (!requireNamespace("BEDMatrix", quietly = TRUE))
      stop("reading a .bed file needs the 'BEDMatrix' package", call. = FALSE)
    path <- sub("\\.bed$", "", genotype)
    G <- BEDMatrix::BEDMatrix(path = path, simple_names = TRUE)
    bim_file <- paste0(path, ".bim")
    bim <- if (file.exists(bim_file)) {
      b <- utils::read.table(bim_file, header = FALSE, colClasses = c("character", "character",
                             "numeric", "numeric", "character", "character"))
      ## BEDMatrix counts the first allele listed in the .bim (column 5)
      data.frame(variant = b[[2]], chr = b[[1]], pos = b[[4]], effect_allele = b[[5]],
                 other_allele = b[[6]], stringsAsFactors = FALSE)
    }
    return(list(G = G, info = bim))
  }
  if (is.null(rownames(genotype)) || is.null(colnames(genotype)))
    stop("genotype must have sample IDs as rownames and variant IDs as colnames", call. = FALSE)
  list(G = genotype, info = NULL)
}

#' Genome-Wide SynPALM Score Tests
#'
#' Runs the score tests of every method in \code{null} over a set of variants,
#' in chunks.
#'
#' Genotypes are used as additive allele counts. As in the UK Biobank
#' analysis, missing genotypes are set to 0 and variants that are constant
#' among the observed independent individuals are perturbed with
#' \code{\link{fix_constant_columns}} so that the test is defined.
#'
#' @param null Output of \code{\link{synpalm_null}}.
#' @param genotype Path to a PLINK \code{.bed} file (read with
#'   \pkg{BEDMatrix}; chromosome, position and alleles are taken from the
#'   \code{.bim}), or any matrix-like object supporting \code{G[i, j]} with
#'   sample IDs as rownames and variant IDs as colnames.
#' @param variants Variant IDs (or column indices) to test; default all. Use
#'   this to split a genome-wide scan into parallel jobs.
#' @param chunk_size Variants read and tested at a time.
#' @param verbose Print progress.
#'
#' @return A data frame with one row per variant: \code{variant}, and when a
#'   \code{.bim} is available \code{chr}, \code{pos}, \code{effect_allele}
#'   (the counted allele), \code{other_allele}; \code{af} (frequency of the
#'   counted allele in the analysis set, among non-missing calls) and
#'   \code{n_missing}; missing genotypes are mean-imputed with \code{2 * af}
#'   before testing; then for each
#'   method \code{<method>_beta}, \code{<method>_se}, \code{<method>_p} and
#'   \code{<method>_log10p} (\eqn{-\log_{10} p}).
#' @export
synpalm_scan <- function(null, genotype, variants = NULL, chunk_size = 200L, verbose = TRUE) {
  stopifnot(inherits(null, "synpalm_null"))
  src <- .open_genotype(genotype)
  G   <- src$G

  matched <- match(null$ids, as.character(rownames(G)))
  if (anyNA(matched)) stop(sum(is.na(matched)), " analysis individuals are missing from the genotype data",
                           call. = FALSE)
  if (is.null(variants)) variants <- colnames(G)
  if (is.numeric(variants)) variants <- colnames(G)[variants]
  miss_v <- setdiff(variants, colnames(G))
  if (length(miss_v)) stop(length(miss_v), " variants not found in the genotype data", call. = FALSE)

  fix_idx <- intersect(null$obs_protein_index, null$independent_indices)
  chunks  <- split(variants, ceiling(seq_along(variants) / chunk_size))
  out <- vector("list", length(chunks))

  for (ci in seq_along(chunks)) {
    v  <- chunks[[ci]]
    Gm <- as.matrix(G[matched, v, drop = FALSE])
    colnames(Gm) <- v
    n_missing <- colSums(is.na(Gm))
    af <- colMeans(Gm, na.rm = TRUE) / 2
    ## mean imputation: a missing genotype gets the variant's mean dosage 2 * af
    na_at <- which(is.na(Gm), arr.ind = TRUE)
    if (nrow(na_at)) Gm[na_at] <- ifelse(is.nan(af), 0, 2 * af)[na_at[, 2]]
    Gm <- fix_constant_columns(Gm, fix_idx)

    res <- data.frame(variant = v, af = unname(af), n_missing = unname(n_missing),
                      stringsAsFactors = FALSE)
    for (meth in null$methods) {
      r <- .score_fun[[meth]](Gm, null$pars[[meth]], null$independent_indices)
      b  <- as.numeric(r[[paste0("hat_beta_", meth)]])
      vb <- as.numeric(r[[paste0("var_hat_beta_", meth)]])
      lp <- as.numeric(r[[paste0("negative_log10_pval_", meth)]])
      res[[paste0(meth, "_beta")]]   <- b
      res[[paste0(meth, "_se")]]     <- sqrt(vb)
      res[[paste0(meth, "_p")]]      <- 10^(-lp)
      res[[paste0(meth, "_log10p")]] <- lp
    }
    out[[ci]] <- res
    if (verbose && (ci %% 10 == 0 || ci == length(chunks)))
      message(sprintf("scanned %d / %d variants", min(ci * chunk_size, length(variants)), length(variants)))
  }
  out <- do.call(rbind, out)
  if (!is.null(src$info)) {
    info <- src$info[match(out$variant, src$info$variant), -1, drop = FALSE]
    out  <- cbind(out["variant"], info, out[setdiff(names(out), "variant")])
  }
  rownames(out) <- NULL
  out
}

## ===========================================================================
## One call
## ===========================================================================

#' SynPALM GWAS of One Protein, End to End
#'
#' Builds relatedness-aware folds, predicts the protein with cross-fitted
#' random forests, fits the null model and scans the genotypes.
#'
#' The analysis set is the individuals present in \code{grm},
#' \code{covariates}, \code{rf_features} (and \code{rf_binary} if given), in
#' the row order of \code{grm}. Every one of them must also be present in the
#' genotype data. \code{protein} does not define the analysis set: it may hold
#' only the measured individuals (as a biobank protein table usually does), and
#' everyone in the analysis set who is absent from it, or \code{NA} in it,
#' counts as unmeasured and contributes through the synthetic phenotype.
#'
#' @param protein Named numeric vector of the protein, named by sample ID.
#'   Individuals that are \code{NA} or absent are unmeasured.
#' @param covariates Data frame of LMM adjustment covariates, rownames = sample IDs.
#' @param rf_features Data frame of random-forest predictors, rownames = sample IDs.
#' @param grm Sparse GRM with sample IDs as dimnames.
#' @param genotype PLINK \code{.bed} path or matrix-like genotype object; see
#'   \code{\link{synpalm_scan}}.
#' @param rf_binary,rf_fixed,screen,n_top,binary_sig,binary_top,num.trees,num.threads
#'   Passed to \code{\link{synpalm_predict}}.
#' @param K Number of cross-fitting folds.
#' @param methods,vc_constraint Passed to \code{\link{synpalm_null}}.
#' @param variants,chunk_size Passed to \code{\link{synpalm_scan}}.
#' @param folds Optional precomputed \code{\link{synpalm_folds}} result; must
#'   match the analysis set.
#' @param verbose Print progress.
#'
#' @return A list of class \code{"synpalm_gwas"} with \code{results} (the data
#'   frame from \code{\link{synpalm_scan}} plus a column \code{rho_oof}, the
#'   cross-fitted correlation between synthetic and measured protein),
#'   \code{accuracy} (per-fold and overall correlation, see
#'   \code{\link{synpalm_predict}}), \code{prediction}, \code{null} and
#'   \code{folds}.
#' @export
synpalm_gwas <- function(protein, covariates, rf_features, grm, genotype,
                         rf_binary = NULL, rf_fixed = NULL, K = 5L,
                         methods = c("SynSurrG", "ObsG"), vc_constraint = "none",
                         variants = NULL, chunk_size = 200L, folds = NULL,
                         screen = TRUE, n_top = 100L, binary_sig = 5, binary_top = c(5L, 10L),
                         num.trees = 300L, num.threads = 1L, verbose = TRUE) {
  grm <- .check_grm(grm)
  keep <- Reduce(intersect, Filter(Negate(is.null), list(
    .ids_of(covariates, "covariates"), .ids_of(rf_features, "rf_features"),
    if (!is.null(rf_binary)) .ids_of(rf_binary, "rf_binary"))))
  ids <- rownames(grm)[rownames(grm) %in% keep]
  if (!length(ids)) stop("no individuals in common between grm and the inputs", call. = FALSE)
  grm <- grm[ids, ids]

  ## the protein over the whole analysis set: absent = unmeasured
  mp <- match(ids, .ids_of(protein, "protein"))
  protein <- setNames(as.numeric(protein)[mp], ids)
  if (verbose) message(sprintf("analysis set: %d individuals, %d with a measured protein",
                               length(ids), sum(!is.na(protein))))

  if (is.null(folds)) {
    folds <- synpalm_folds(grm, K = K)
  } else if (!identical(folds$ids, ids)) {
    stop("folds$ids does not match the analysis set", call. = FALSE)
  }

  pred <- synpalm_predict(protein, rf_features, folds, rf_binary = rf_binary, rf_fixed = rf_fixed,
                          screen = screen, n_top = n_top, binary_sig = binary_sig,
                          binary_top = binary_top, num.trees = num.trees,
                          num.threads = num.threads, verbose = verbose)
  if (verbose) message("fitting the null model")
  null <- synpalm_null(pred$protein, pred$protein_hat, covariates, grm, folds,
                       methods = methods, vc_constraint = vc_constraint)
  res <- synpalm_scan(null, genotype, variants = variants, chunk_size = chunk_size, verbose = verbose)
  res$rho_oof <- pred$rho_oof

  structure(list(results = res, accuracy = pred$accuracy, prediction = pred, null = null,
                 folds = folds),
            class = "synpalm_gwas")
}

#' @method print synpalm_prediction
#' @export
print.synpalm_prediction <- function(x, digits = 4, ...) {
  K <- nrow(x$accuracy) - 1L
  cat(sprintf("SynPALM synthetic phenotype: %d individuals, %d with a measured protein\n",
              length(x$ids), x$n_obs))
  cat(sprintf("%d-fold cross-fitted random forest; correlation with the measured protein:\n", K))
  print(x$accuracy, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @method print synpalm_gwas
#' @export
print.synpalm_gwas <- function(x, digits = 4, ...) {
  print(x$prediction, digits = digits)
  cat(sprintf("\nGWAS: %d variants, methods %s; see $results\n", nrow(x$results),
              paste(x$null$methods, collapse = ", ")))
  invisible(x)
}

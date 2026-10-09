options(synsurrg.verbose = FALSE)

test_that("folds keep relatedness clusters intact and leave the RNG alone", {
  d <- sim_cohort(n_fam = 300)
  set.seed(99); before <- runif(1); set.seed(99)
  f <- synpalm_folds(d$grm, K = 5)
  expect_equal(runif(1), before)
  expect_equal(f$ids, rownames(d$grm))
  expect_true(all(vapply(f$blocks, function(ix) length(unique(f$fold[ix])) == 1, logical(1))))
  expect_lt(diff(range(table(f$fold))), 0.05 * length(f$fold))
  expect_identical(f, synpalm_folds(d$grm, K = 5))
})

test_that("pipeline steps reproduce the low-level functions exactly", {
  d <- sim_cohort()
  f <- synpalm_folds(d$grm)
  pr <- synpalm_predict(d$protein, d$rf_features, f, rf_binary = d$rf_binary,
                        rf_fixed = "age", num.trees = 50, verbose = FALSE)
  expect_false(anyNA(pr$protein_hat))
  expect_gt(pr$rho_oof, 0.3)

  nl <- synpalm_null(pr$protein, pr$protein_hat, d$covariates, d$grm, f,
                     methods = c("SynSurrG", "ObsG", "SynSurr", "Obs"))
  res <- synpalm_scan(nl, d$G, verbose = FALSE)

  ## the same model fitted by hand
  X <- model.matrix(~ ., d$covariates)[, -1]
  colnames(X) <- make.names(colnames(X))
  n <- length(f$ids)
  r <- rank(pr$protein_hat)
  mydf <- list(X_all = as.data.frame(X),
               S = qnorm((r - 0.375) / (n - 2 * 0.375 + 1)),
               Y_obs = INT(data.frame(protein = pr$protein, row = seq_len(n)), "protein")$int,
               GRM = d$grm)
  old <- options(synsurrg.on_nonpsd = "none")
  p1 <- SynSurrG_ablation_estimate(mydf)
  options(old)
  Gm <- fix_constant_columns(d$G, intersect(which(!is.na(pr$protein)), f$independent_indices))
  t1 <- score_test_SynSurrG_multiply(Gm, p1)
  expect_equal(res$SynSurrG_log10p, as.numeric(t1$negative_log10_pval_SynSurrG))
  expect_equal(res$SynSurrG_beta, as.numeric(t1$hat_beta_SynSurrG))
  expect_equal(res$SynSurrG_se, sqrt(as.numeric(t1$var_hat_beta_SynSurrG)))
  expect_true(all(c("ObsG_p", "SynSurr_p", "Obs_p") %in% names(res)))
  expect_equal(nl$pars$SynSurrG$vc$mode, "none")
})

test_that("synpalm_gwas aligns inputs by sample ID, not by position", {
  d <- sim_cohort(n_fam = 400, n_snp = 10)
  a <- synpalm_gwas(d$protein, d$covariates, d$rf_features, d$grm, d$G,
                    rf_binary = d$rf_binary, num.trees = 50, verbose = FALSE)
  o <- sample(length(d$ids))
  b <- synpalm_gwas(d$protein[o], d$covariates[o, ], d$rf_features[rev(o), ], d$grm,
                    d$G[o, ], rf_binary = d$rf_binary[o, ], num.trees = 50, verbose = FALSE)
  expect_equal(a$results, b$results)
  expect_gt(a$results$SynSurrG_log10p[1], a$results$ObsG_log10p[1])
})

test_that("vc_constraint switches between raw and projected estimates", {
  d <- sim_cohort(n_fam = 300)
  f <- synpalm_folds(d$grm)
  pr <- synpalm_predict(d$protein, d$rf_features, f, num.trees = 30, verbose = FALSE)
  a <- synpalm_null(pr$protein, pr$protein_hat, d$covariates, d$grm, f, methods = "SynSurrG")
  b <- synpalm_null(pr$protein, pr$protein_hat, d$covariates, d$grm, f, methods = "SynSurrG",
                    vc_constraint = "project")
  expect_equal(a$pars$SynSurrG$vc$raw, b$pars$SynSurrG$vc$raw)
  expect_equal(a$pars$SynSurrG$vc$used[-4], a$pars$SynSurrG$vc$raw[-4])
  expect_equal(getOption("synsurrg.on_nonpsd"), "none")
})

test_that("PLINK .bed input matches the matrix input and carries .bim info", {
  skip_if_not_installed("BEDMatrix")
  d <- sim_cohort(n_fam = 300, n_snp = 8)
  G <- d$G; G[3, 2] <- NA
  prefix <- file.path(tempdir(), "synpalm_test")
  write_plink(G, prefix)
  f  <- synpalm_folds(d$grm)
  pr <- synpalm_predict(d$protein, d$rf_features, f, num.trees = 30, verbose = FALSE)
  nl <- synpalm_null(pr$protein, pr$protein_hat, d$covariates, d$grm, f)
  a <- synpalm_scan(nl, G, verbose = FALSE)
  b <- synpalm_scan(nl, paste0(prefix, ".bed"), verbose = FALSE)
  expect_equal(b[names(a)], a)
  expect_equal(b$effect_allele, rep("A", 8))
  expect_equal(b$pos, (1:8) * 100)
  expect_equal(a$n_missing[2], 1)
  Gi <- G; Gi[3, 2] <- mean(G[, 2], na.rm = TRUE)
  ai <- synpalm_scan(nl, Gi, verbose = FALSE)
  expect_equal(a[setdiff(names(a), "n_missing")], ai[setdiff(names(ai), "n_missing")])
  sub <- synpalm_scan(nl, paste0(prefix, ".bed"), variants = c("rs7", "rs2"), verbose = FALSE)
  expect_equal(sub$SynSurrG_log10p, a$SynSurrG_log10p[c(7, 2)])
})

test_that("folds reorder a GRM whose relatives are scattered", {
  d <- sim_cohort(n_fam = 300)
  o <- sample(nrow(d$grm))
  scrambled <- d$grm[o, o]
  f <- synpalm_folds(scrambled)
  g <- synpalm_folds(d$grm)
  expect_equal(length(f$blocks), 300)
  expect_setequal(f$ids, g$ids)
  expect_lt(max(lengths(f$blocks)), 5)
  ## an already contiguous order is kept as is
  expect_identical(g$ids, rownames(d$grm))
  ## the whole pipeline still runs and finds the causal variant
  r <- synpalm_gwas(d$protein, d$covariates, d$rf_features, scrambled, d$G[, 1:5],
                    num.trees = 50, verbose = FALSE)
  expect_gt(r$results$SynSurrG_log10p[1], 2)
})

test_that("5-fold cross-fit: a fold's predictions never use its own phenotypes", {
  d <- sim_cohort(n_fam = 300)
  f <- synpalm_folds(d$grm)
  expect_equal(f$K, 5L)
  pr <- synpalm_predict(d$protein, d$rf_features, f, rf_binary = d$rf_binary,
                        num.trees = 50, verbose = FALSE)

  ## scramble the measured proteins of fold 1: fold 1 predictions must not move
  y2 <- d$protein
  in1 <- f$ids[f$fold == 1]
  y2[in1] <- rev(y2[in1])
  pr2 <- synpalm_predict(y2, d$rf_features, f, rf_binary = d$rf_binary,
                         num.trees = 50, verbose = FALSE)
  expect_identical(pr2$protein_hat[f$fold == 1], pr$protein_hat[f$fold == 1])
  expect_false(identical(pr2$protein_hat[f$fold != 1], pr$protein_hat[f$fold != 1]))

  ## the accuracy table is the correlation computed directly
  acc <- pr$accuracy
  expect_equal(acc$fold, c(as.character(1:5), "all"))
  obs <- !is.na(pr$protein)
  expect_equal(acc$rho[6], cor(pr$protein[obs], pr$protein_hat[obs]))
  expect_equal(acc$rho[6], pr$rho_oof)
  for (k in 1:5) {
    ii <- obs & f$fold == k
    expect_equal(acc$rho[k], cor(pr$protein[ii], pr$protein_hat[ii]))
    expect_equal(acc$n_labelled[k], sum(ii))
    expect_equal(acc$n_train[k], sum(obs & f$fold != k))
  }
  expect_equal(sum(acc$n_predicted[1:5]), length(f$ids))
  expect_equal(acc$r2, acc$rho^2)

  g <- suppressWarnings(synpalm_gwas(d$protein, d$covariates, d$rf_features, d$grm, d$G[, 1:3],
                                     rf_binary = d$rf_binary, num.trees = 50, verbose = FALSE))
  expect_true(all(g$results$rho_oof == g$accuracy$rho[6]))
  expect_output(print(g), "5-fold cross-fitted")
})

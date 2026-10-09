## Small simulated cohort: families of 1-4 relatives, a protein measured in
## `obs_frac` of individuals, a synthetic-phenotype feature set, covariates
## and genotypes. Variant 1 is causal.
sim_cohort <- function(n_fam = 600, obs_frac = 0.3, n_snp = 40, beta_g = 0.3, seed = 1) {
  set.seed(seed)
  sz  <- sample(1:4, n_fam, replace = TRUE, prob = c(.5, .2, .2, .1))
  n   <- sum(sz)
  ids <- sprintf("id%05d", seq_len(n))
  fam <- rep(seq_len(n_fam), sz)
  ## GRM: 0.5 within families, 1 on the diagonal
  pairs <- do.call(rbind, lapply(split(seq_len(n), fam), function(ix) expand.grid(i = ix, j = ix)))
  grm <- Matrix::sparseMatrix(i = pairs$i, j = pairs$j,
                              x = ifelse(pairs$i == pairs$j, 1, 0.5), dims = c(n, n),
                              dimnames = list(ids, ids))
  L <- Matrix::t(Matrix::chol(grm))
  u <- as.vector(L %*% rnorm(n)) * sqrt(0.3)

  covariates <- data.frame(age = rnorm(n), sex = factor(sample(c("F", "M"), n, TRUE)),
                           PC1 = rnorm(n), PC2 = rnorm(n), row.names = ids)
  G <- vapply(seq_len(n_snp), function(j) rbinom(n, 2, runif(1, 0.1, 0.4)), numeric(n))
  dimnames(G) <- list(ids, sprintf("rs%d", seq_len(n_snp)))

  y <- 0.2 * covariates$age + beta_g * scale(G[, 1])[, 1] + u + rnorm(n, sd = sqrt(0.6))
  rf_features <- data.frame(f1 = y + rnorm(n), f2 = 0.5 * y + rnorm(n), f3 = rnorm(n),
                            f4 = rnorm(n), age = covariates$age, row.names = ids)
  rf_binary <- data.frame(d1 = rbinom(n, 1, plogis(y)), d2 = rbinom(n, 1, 0.3), row.names = ids)

  protein <- setNames(rep(NA_real_, n), ids)
  obs <- sample.int(n, floor(obs_frac * n))
  protein[obs] <- y[obs]
  list(protein = protein, covariates = covariates, rf_features = rf_features,
       rf_binary = rf_binary, grm = grm, G = G, ids = ids)
}

## Minimal PLINK .bed/.bim/.fam writer (SNP-major). Genotypes are counts of
## the first .bim allele, which is what BEDMatrix returns.
write_plink <- function(G, prefix) {
  n <- nrow(G); p <- ncol(G)
  code <- function(g) ifelse(is.na(g), 1L, ifelse(g == 2, 0L, ifelse(g == 1, 2L, 3L)))
  nb <- ceiling(n / 4)
  bytes <- unlist(lapply(seq_len(p), function(j) {
    c4 <- c(code(G[, j]), rep(0L, nb * 4 - n))
    m <- matrix(c4, nrow = 4)
    as.raw(m[1, ] + 4L * m[2, ] + 16L * m[3, ] + 64L * m[4, ])
  }))
  writeBin(c(as.raw(c(0x6c, 0x1b, 0x01)), bytes), paste0(prefix, ".bed"))
  write.table(data.frame(rownames(G), rownames(G), 0, 0, 0, -9),
              paste0(prefix, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(data.frame(1, colnames(G), 0, seq_len(p) * 100, "A", "G"),
              paste0(prefix, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(prefix)
}

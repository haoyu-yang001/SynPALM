## General data-generating model for the new simulations of the revision
## (S2 alpha_G-only SNPs, S3/S6 genetic vs residual sharing, S4 proportion of
## related individuals, S5 MAR missingness). Sourced after code/funcs_v2.R.
##
## Cohort: n_all individuals, a fraction prop_related of them in pairs with
## kinship r (individuals 2i-1, 2i), the rest unrelated; GRM K block diagonal.
## For a SNP g and covariate x (standard normal, "age"):
##   Y = beta_G g + beta_X x + u_T + e_T,   S = alpha_G g + alpha_X x + u_S + e_S,
##   (u_T, u_S) ~ N(0, Gmat (x) K),  (e_T, e_S) ~ N(0, Emat (x) I),
## Gmat = [tauT^2, tauTS; tauTS, tauS^2], Emat = [sigmaT^2, sigmaTS; sigmaTS, sigmaS^2]
## (tau / sigma are standard deviations, as in funcs_v2.R). Within a pair the
## random effects are drawn as sqrt(r) f + sqrt(1 - r) w_i with f, w_i ~ N(0, Gmat),
## which gives exactly Cov = r Gmat. Y is then set to NA for the unlabelled.

## covariance matrices from genetic and residual correlations
vc_from_rho <- function(rho_g, rho_e, tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5) {
  stopifnot(abs(rho_g) < 1, abs(rho_e) < 1)
  Gm <- matrix(c(tauT^2, rho_g * tauT * tauS, rho_g * tauT * tauS, tauS^2), 2)
  Em <- matrix(c(sigmaT^2, rho_e * sigmaT * sigmaS, rho_e * sigmaT * sigmaS, sigmaS^2), 2)
  list(Gm = Gm, Em = Em, rho = (Gm[1, 2] + Em[1, 2]) / sqrt((Gm[1, 1] + Em[1, 1]) * (Gm[2, 2] + Em[2, 2])),
       rho_g = rho_g, rho_e = rho_e)
}
## the submitted 3/5 : 2/5 split (with the rho = 0.75 correction)
vc_from_total <- function(rho, tauT = 0.7, tauS = 0.4, sigmaT = 0.7, sigmaS = 0.5) {
  cv <- sim_cov_TS(rho, tauT, tauS, sigmaT, sigmaS)
  vc_from_rho(cv[["rho_g"]], cv[["rho_e"]], tauT, tauS, sigmaT, sigmaS)
}

make_cohort <- function(n_all, prop_related = 1, r = 0.5) {
  n_pair <- 2L * round(n_all * prop_related / 2)
  fam <- c(rep(seq_len(n_pair / 2), each = 2L), n_pair / 2 + seq_len(n_all - n_pair))
  i <- seq_len(n_pair)
  odd <- i[i %% 2 == 1]
  GRM <- sparseMatrix(i = c(seq_len(n_all), odd, odd + 1L), j = c(seq_len(n_all), odd + 1L, odd),
                      x = c(rep(1, n_all), rep(r, 2 * length(odd))), dims = c(n_all, n_all))
  list(n_all = n_all, n_pair = n_pair, r = r, fam = fam, GRM = GRM,
       pair1 = odd, pair2 = odd + 1L, single = if (n_pair < n_all) (n_pair + 1L):n_all else integer(0))
}

## random effects + residuals (n_all x 2: T, S)
draw_noise <- function(coh, vc) {
  n <- coh$n_all; cG <- chol(vc$Gm); cE <- chol(vc$Em)
  u <- matrix(rnorm(2 * n), n) %*% cG                     # w_i (and singletons)
  if (coh$n_pair) {
    f <- matrix(rnorm(coh$n_pair), coh$n_pair / 2) %*% cG # shared family part
    k <- c(coh$pair1, coh$pair2)
    u[k, ] <- sqrt(1 - coh$r) * u[k, ] + sqrt(coh$r) * rbind(f, f)
  }
  u + matrix(rnorm(2 * n), n) %*% cE
}

## gene drop consistent with the cohort: pairs are siblings (r = 0.5),
## half-siblings (0.25) or first cousins (0.125); singletons are founders
gene_drop_cohort <- function(coh, m, maf) {
  n <- coh$n_all; np <- coh$n_pair / 2
  G <- matrix(rbinom(n * m, 2L, maf), n, m)
  if (np) {
    founder <- function() matrix(rbinom(np * m, 2L, maf), np, m)
    gam <- function(g) matrix(rbinom(length(g), 1L, g / 2), nrow(g))
    gf <- founder(); gm <- founder()
    pr <- if (coh$r == 0.5) list(gam(gf) + gam(gm), gam(gf) + gam(gm))
          else if (coh$r == 0.25) list(gam(gf) + gam(gm), gam(gf) + gam(founder()))
          else { f <- gam(gf) + gam(gm); u <- gam(gf) + gam(gm)
                 list(gam(f) + gam(founder()), gam(u) + gam(founder())) }
    G[coh$pair1, ] <- pr[[1]]; G[coh$pair2, ] <- pr[[2]]
  }
  G
}

## labelled set
##   "pairs"  : whole families labelled together, first families first (submitted design)
##   "random" : individuals labelled completely at random
##   "mar_x"  : P(labelled) = plogis(a + b x), a solved so the labelled fraction is 1 - miss
##   "mar_s"  : the same with S (the synthetic phenotype) in place of x
label_set <- function(coh, miss, type = "pairs", x = NULL, s = NULL, b = 1) {
  n <- coh$n_all; n_lab <- round(n * (1 - miss))
  if (miss == 0) return(seq_len(n))
  switch(type,
    pairs  = { o <- unique(coh$fam); keep <- o[cumsum(tabulate(coh$fam)[o]) <= n_lab]
               which(coh$fam %in% keep) },
    random = sort(sample.int(n, n_lab)),
    mar_x  = , mar_s = {
      z <- if (type == "mar_x") x else s; z <- (z - mean(z)) / sd(z)
      a <- uniroot(function(a) mean(plogis(a + b * z)) - (1 - miss), c(-20, 20))$root
      which(runif(n) < plogis(a + b * z)) })
}

## mydf in the format of the *_typeIerror_step1 functions of funcs_v2.R
make_mydf <- function(coh, Y, S, x, lab) {
  Y_obs <- rep(NA_real_, coh$n_all); Y_obs[lab] <- Y[lab]
  list(X_all = x, S = S, Y = Y, Y_obs = Y_obs, GRM = coh$GRM)
}

## one individual per family (for SynSurr / Obs, which assume independence)
one_per_family <- function(coh, seed = 25324) withr::with_seed(seed,
  as.integer(tapply(seq_len(coh$n_all), coh$fam, function(i) i[sample.int(length(i), 1L)])))

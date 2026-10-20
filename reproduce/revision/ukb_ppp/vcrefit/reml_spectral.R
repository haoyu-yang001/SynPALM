## Exact REML for the Syn-PALM variance components with a sparse (block-diagonal) GRM.
## The sparse GRM K is block diagonal over the relatedness components, so its eigen-
## decomposition K = U diag(lambda) U' is computed block by block (largest UKB block
## ~4,700). After rotation by U' the covariance of one trait is diag(tau lambda + sigma),
## and of two traits observed on the same individuals block-diag(lambda_i G + E), so the
## REML log-likelihood is evaluated exactly in O(n p^2) per evaluation.
##   reml_uni(): (tau, sigma) of one trait; profiled over the scale, 1-D in h = tau/(tau+sigma),
##               h in [0, h_max) (tau >= 0, the REML parameter space)
##   reml_biv(): 2x2 G and E of two traits on the same individuals, G = L L', E = M M'
##               (positive semi-definite by construction), BFGS from the univariate fits
## Used by code/post_cf3/vcrefit/vcrefit.R. GMMAT::glmmkin (the REML of the VC
## simulations) is used as a check of reml_uni(); its multi-phenotype fit returned zero
## cross-covariances (boundary refit) for the UKB data and is not used.
spec_decomp <- function(K, sub, memb) {
  ## K: n x n sparse GRM; sub: indices of the individuals used; memb: component id of each
  ## of the n individuals (K[i, j] != 0 only within a component)
  Ks <- as(as(K[sub, sub], "generalMatrix"), "TsparseMatrix")
  ns <- length(sub); ms <- memb[sub]
  lam <- numeric(ns)
  grp <- split(seq_len(ns), ms)
  sz  <- lengths(grp)
  dg  <- diag(Ks)
  single <- unlist(grp[sz == 1L], use.names = FALSE)
  lam[single] <- dg[single]
  multi <- grp[sz > 1L]
  ii <- Ks@i + 1L; jj <- Ks@j + 1L; xx <- Ks@x
  keep <- ms[ii] %in% names(multi)
  stopifnot(all(ms[ii] == ms[jj]))                 # block diagonal
  trip <- split(seq_along(ii)[keep], ms[ii][keep])
  loc <- integer(ns)
  U <- vector("list", length(multi))
  for (b in seq_along(multi)) {
    g <- multi[[b]]; loc[g] <- seq_along(g)
    t <- trip[[names(multi)[b]]]
    M <- matrix(0, length(g), length(g)); M[cbind(loc[ii[t]], loc[jj[t]])] <- xx[t]
    e <- eigen(M, symmetric = TRUE)
    lam[g] <- e$values; U[[b]] <- e$vectors
  }
  list(lam = lam, groups = multi, U = U, n = ns)
}
spec_rotate <- function(S, M) {                    # U' M (rows = individuals of `sub`)
  M <- as.matrix(M); out <- M
  for (b in seq_along(S$groups)) { g <- S$groups[[b]]; out[g, ] <- crossprod(S$U[[b]], M[g, , drop = FALSE]) }
  out
}
reml_uni <- function(S, y, X) {
  ys <- as.numeric(spec_rotate(S, y)); Xs <- spec_rotate(S, X)
  n <- length(ys); p <- qr(X)$rank; lam <- S$lam
  nll <- function(h) {
    d <- h * lam + (1 - h); w <- 1 / d
    A <- crossprod(Xs, Xs * w); b <- crossprod(Xs, ys * w)
    ch <- chol(A); beta <- backsolve(ch, forwardsolve(t(ch), b))
    r <- ys - Xs %*% beta; Q <- sum(w * r^2)
    0.5 * ((n - p) * log(Q / (n - p)) + sum(log(d)) + 2 * sum(log(diag(ch))))
  }
  hmax <- if (min(lam) < 1) min(0.999, (1 - 1e-6) / (1 - min(lam))) else 0.999
  o <- optimize(nll, c(0, hmax), tol = 1e-8)
  h <- if (nll(0) <= o$objective) 0 else o$minimum
  d <- h * lam + (1 - h); w <- 1 / d
  A <- crossprod(Xs, Xs * w); beta <- solve(A, crossprod(Xs, ys * w))
  s <- sum(w * (ys - Xs %*% beta)^2) / (n - p)
  list(tau = h * s, sigma = (1 - h) * s, h = h, nll = min(o$objective, nll(0)), boundary = h == 0,
       rotated = list(ys = ys, Xs = Xs))
}
reml_biv <- function(S, y, s, X, start = NULL, rot_y = NULL, rot_s = NULL, rot_X = NULL, psd = TRUE) {
  ## psd = TRUE : G = L L', E = M M' (REML on the parameter space)
  ## psd = FALSE: the six entries free, only every lambda_i G + E positive definite
  ##              (where the unconstrained REML optimum lies; start = list(G, E) required)
  ys <- if (is.null(rot_y)) as.numeric(spec_rotate(S, y)) else rot_y
  ss <- if (is.null(rot_s)) as.numeric(spec_rotate(S, s)) else rot_s
  Xs <- if (is.null(rot_X)) spec_rotate(S, X) else rot_X
  lam <- S$lam; p <- ncol(Xs)
  mats <- function(th) {
    if (!psd) return(list(G = matrix(th[c(1, 2, 2, 3)], 2), E = matrix(th[c(4, 5, 5, 6)], 2)))
    L <- matrix(c(th[1], th[2], 0, th[3]), 2); M <- matrix(c(th[4], th[5], 0, th[6]), 2)
    list(G = L %*% t(L), E = M %*% t(M))
  }
  nll <- function(th) {
    m <- mats(th); G <- m$G; E <- m$E
    a <- lam * G[1, 1] + E[1, 1]; b <- lam * G[1, 2] + E[1, 2]; c <- lam * G[2, 2] + E[2, 2]
    dt <- a * c - b^2
    if (any(!is.finite(dt)) || any(dt <= 0) || any(a <= 0)) return(1e300)
    w11 <- c / dt; w12 <- -b / dt; w22 <- a / dt
    A <- rbind(cbind(crossprod(Xs, Xs * w11), crossprod(Xs, Xs * w12)),
               cbind(crossprod(Xs, Xs * w12), crossprod(Xs, Xs * w22)))
    v <- c(crossprod(Xs, w11 * ys + w12 * ss), crossprod(Xs, w12 * ys + w22 * ss))
    ch <- tryCatch(chol(A), error = function(e) NULL); if (is.null(ch)) return(1e300)
    beta <- backsolve(ch, forwardsolve(t(ch), v))
    ry <- ys - Xs %*% beta[1:p]; rs <- ss - Xs %*% beta[p + 1:p]
    Q <- sum(w11 * ry^2 + 2 * w12 * ry * rs + w22 * rs^2)
    0.5 * (sum(log(dt)) + 2 * sum(log(diag(ch))) + Q)
  }
  th0 <- function(G, E) {                             # Cholesky of start values (made PD)
    f <- function(A) { e <- eigen(A, symmetric = TRUE); e$values <- pmax(e$values, 1e-4 * max(e$values, 1e-4))
      t(chol(e$vectors %*% diag(e$values) %*% t(e$vectors))) }
    L <- f(G); M <- f(E); c(L[1, 1], L[2, 1], L[2, 2], M[1, 1], M[2, 1], M[2, 2])
  }
  starts <- if (psd) list(th0(start$G, start$E), th0(diag(diag(start$G) + 1e-3), start$E))   # 2nd start: rho_g = 0
            else list(c(start$G[c(1, 2, 4)], start$E[c(1, 2, 4)]))
  fits <- lapply(starts, function(t0) {
    f <- optim(t0, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12, ndeps = rep(1e-5, 6)))
    ## polish with Nelder-Mead (robust where the PD boundary of the free parametrisation is near)
    g <- optim(f$par, nll, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-14))
    if (g$value < f$value) { g$convergence <- f$convergence; g } else f })
  best <- fits[[which.min(vapply(fits, `[[`, 0, "value"))]]
  m <- mats(best$par)
  list(G = m$G, E = m$E, nll = best$value, convergence = best$convergence,
       nll_starts = vapply(fits, `[[`, 0, "value"))
}

## =====================================================================
## funcs_fastscan.R -- vectorised Syn-PALM / SynSurr score test.
##
## score_test_SynSurrG_multiply() in funcs_v2.R loops over SNPs and, for
## every SNP, solves the (1+p) x (1+p) system A_mat x = [Ah_d; Atb] with an
## n_obs-column right-hand side. A_mat differs between SNPs only in its first
## row and column, so the solve has a closed form through the Schur
## complement s_d = A11_d - A12_d A22^{-1} A12_d'. Writing it out removes the
## loop; every remaining operation is a matrix product over the whole chunk.
##
##   Q   = A22^{-1}                      (p x p, shared by all SNPs)
##   alpha1_d = (h_d - c_d Q w) / s_d,  alpha2_d = Q w - Q c_d' alpha1_d
##   S_d  = Stt1_d . e0 + alpha1_d hh1_d1 + hh1_d,-1 . alpha2_d
##   V2_d = VV2_1_d . VV1_d + (hh1_d1 - k_d)(u_d - c_d Q T_d)/s_d + hh1_d,-1 Q T_d
## with c_d = A12[d,], h_d = G_d' Sigma22^{-1} hatY, w = Att' hatY,
## e0 = Y - X_obs (B1 - B2) - B2_1, u_d = Ah_d . VV1_d, T = Atb VV1',
## k_d = hh1_d,-1 Q c_d'. Same output as score_test_SynSurrG_multiply().
##
## add_fastscan_cache() precomputes the SNP-invariant pieces once per protein
## (C = Sigma12 Sigma22^{-1} Sigma21, Q, Qw, e0).
## =====================================================================

add_fastscan_cache <- function(pars, Sigma22 = NULL) {
  ## With Sigma22 given, Sigma22^{-1} G is computed by a sparse Cholesky solve.
  ## The explicit inverse of a large connected relatedness component is dense
  ## (nnz ~ 2.8e7 in UKB), its Cholesky factor is not.
  if (!is.null(Sigma22)) {
    pars$ch22 <- Matrix::Cholesky(Matrix::forceSymmetric((Sigma22 + Matrix::t(Sigma22)) / 2),
                                  perm = TRUE, LDL = FALSE)
    pars$inv_Sigma22 <- NULL
  }
  ## Sigma22^{-1} hatY, so that h = G' Sigma22^{-1} hatY needs no n x m product
  pars$u22 <- if (!is.null(pars$ch22)) as.numeric(Matrix::solve(pars$ch22, pars$hatY, system = "A"))
              else as.numeric(pars$inv_Sigma22 %*% pars$hatY)
  pars$Cmat <- as(pars$Sigma12_Sigma22inv %*% Matrix::t(pars$Sigma12), "CsparseMatrix")
  pars$Qinv <- solve(as.matrix(pars$A22))
  pars$Att  <- as.matrix(pars$Att)
  w         <- as.numeric(crossprod(pars$Att, pars$hatY))
  pars$Qw   <- as.numeric(pars$Qinv %*% w)
  pars$e0   <- as.numeric(pars$Y - pars$X_obs %*% (pars$B1 - pars$B2) - pars$B2_1)
  pars
}

score_test_SynSurrG_fast <- function(g_matrix, pars, subset = NULL, label = "SynSurrG") {
  if (!is.null(subset)) g_matrix <- g_matrix[subset, , drop = FALSE]
  G_all <- as.matrix(g_matrix)
  G_obs <- G_all[pars$obs_protein_index, , drop = FALSE]
  Q     <- pars$Qinv

  ## A11_d = g' Sigma22^{-1} g = || L^{-1} P g ||^2 (one triangular solve)
  if (!is.null(pars$ch22)) {
    W   <- as.matrix(Matrix::solve(pars$ch22, Matrix::solve(pars$ch22, G_all, system = "P"), system = "L"))
    A11 <- colSums(W^2); rm(W)
  } else A11 <- colSums(G_all * as.matrix(pars$inv_Sigma22 %*% G_all))
  A12  <- crossprod(G_all, pars$Att)                             # m x p
  h    <- as.numeric(crossprod(G_all, pars$u22))                 # m

  Stt1 <- as.matrix(crossprod(G_obs, pars$V11))                  # m x n_obs
  Stt2 <- as.matrix(crossprod(G_obs, pars$Btt1))                 # m x p
  VV1  <- Stt1 - Stt2 %*% as.matrix(pars$B_mat)                  # m x n_obs
  bt1  <- as.matrix(pars$Sigma12_Sigma22inv %*% G_all)           # n_obs x m
  Ah   <- t(bt1)                    # Sigma12 Sigma22^{-1} G, transposed (m x n_obs)
  bt2  <- as.matrix(pars$bt2)

  hh1_1 <- rowSums(VV1 * t(bt1))                                 # m
  hh1_2 <- VV1 %*% bt2                                           # m x p
  V1    <- rowSums(as.matrix(VV1 %*% pars$Sigma11) * VV1)
  V21   <- -rowSums(as.matrix(VV1 %*% pars$Cmat) * VV1)          # VV2_1 . VV1

  CQ  <- A12 %*% Q                                               # m x p
  s   <- A11 - rowSums(CQ * A12)                                 # Schur complement
  alpha1 <- (h - as.numeric(A12 %*% pars$Qw)) / s
  ## alpha2_d = Qw - Q c_d' alpha1_d  ->  m x p, row d
  alpha2 <- matrix(pars$Qw, nrow(A12), length(pars$Qw), byrow = TRUE) - CQ * alpha1

  S_syn <- as.numeric(Stt1 %*% pars$e0) + alpha1 * hh1_1 + rowSums(hh1_2 * alpha2)

  Tm  <- VV1 %*% t(as.matrix(pars$Atb))                          # m x p  (= (Atb VV1')')
  u   <- rowSums(Ah * VV1)
  HQ  <- hh1_2 %*% Q
  k   <- rowSums(HQ * A12)
  cQT <- rowSums(CQ * Tm)
  hQT <- rowSums(HQ * Tm)
  VU  <- V1 + V21 + (hh1_1 - k) * (u - cQT) / s + hQT

  inv_beta <- rowSums(VV1 * t(G_obs))
  hat_beta <- S_syn / inv_beta
  var_hat  <- VU / inv_beta^2
  sc <- compute_score(S_syn, VU)
  out <- list(sc$score, sc$negative_log10_pval, hat_beta, var_hat)
  names(out) <- paste0(c("T_score_", "negative_log10_pval_", "hat_beta_", "var_hat_beta_"), label)
  out
}

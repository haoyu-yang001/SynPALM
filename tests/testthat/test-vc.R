## the re-solve applies to the "project" mode (vc_constraint = "project");
## the package default "none" keeps the raw estimates as in the submitted version
test_that("residual components are re-solved after the genetic block is constrained", {
  withr::local_options(synsurrg.verbose = FALSE, synsurrg.on_nonpsd = "project")
  d <- 0.9975
  ## negative tau: sigma takes up tau_raw * dbar, so tau * dbar + sigma is kept
  u <- constrain_vc_univariate(-0.16, 1.158, dbar = d)
  expect_equal(u$tau2, 0)
  expect_equal(u$sigma2, 1.158 - 0.16 * d)
  ## valid estimates are untouched
  u <- constrain_vc_univariate(0.3, 0.7, dbar = d)
  expect_equal(c(u$tau2, u$sigma2), c(0.3, 0.7))
  ## bivariate: diagonal moments of T and S are preserved
  b <- constrain_vc_bivariate(-0.05, 0.02, 0.3, 1.05, 0.3, 0.7, dbar_T = d, dbar_S = d)
  expect_true(b$vc$projected_G)
  expect_equal(b$tau_T2 * d + b$sigma_T2, -0.05 * d + 1.05, tolerance = 1e-8)
  expect_equal(b$tau_S2 * d + b$sigma_S2,  0.30 * d + 0.70, tolerance = 1e-8)
  ## PSD-valid genetic block: residual block identical to the old code
  b <- constrain_vc_bivariate(0.4, 0.1, 0.3, 0.6, 0.2, 0.5, dbar_T = d, dbar_S = d)
  expect_equal(unname(unlist(b[1:6])), c(0.4, 0.1, 0.3, 0.6, 0.2, 0.5))
  ## switch back to the 2026-08-22 behaviour
  withr::local_options(synsurrg.rematch_sigma = FALSE)
  expect_equal(constrain_vc_univariate(-0.16, 1.158, dbar = d)$sigma2, 1.158)
})

test_that("no related pairs: the genetic component is set to 0 instead of NaN", {
  withr::local_options(synsurrg.verbose = FALSE)
  expect_equal(SynPALM:::.ss_mom_ratio(1, 0), 0)
  expect_equal(SynPALM:::.ss_mom_ratio(1, 4), 0.25)
})

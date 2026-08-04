## ---------------------------------------------------------------------------
## Replace the 12 placeholder roxygen stubs in R/synpalm_functions.R with real
## documentation. Safe to re-run: it only rewrites the roxygen block directly
## above each named function, and leaves function bodies untouched.
##
## Usage, from the package root:
##   source("tools/fill_placeholder_docs.R")
## ---------------------------------------------------------------------------

TARGET <- "R/synpalm_functions.R"
stopifnot(file.exists(TARGET))

docs <- list(

  reset_temp_dir = c(
    "#' Clear and recreate a scratch directory",
    "#'",
    "#' Deletes \\code{dir} and everything inside it if it exists, then creates it",
    "#' again as an empty directory. Used to give each job a clean temporary",
    "#' workspace before writing intermediate files.",
    "#'",
    "#' @param dir Character path to the directory to reset.",
    "#'",
    "#' @return Called for its side effect on the file system; returns \\code{NULL}",
    "#'   invisibly."),

  ObsG_hatY_typeIerror_step1 = c(
    "#' Null model for the synthetic phenotype among measured individuals",
    "#'",
    "#' Fits the linear mixed null model to the synthetic (predicted) phenotype,",
    "#' restricted to individuals who also have a measured target phenotype, and",
    "#' returns the quantities the corresponding score test reuses across variants.",
    "#' Used to check that the synthetic phenotype alone is well calibrated.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all} (covariates for all",
    "#'   individuals, columns prefixed \"X.\"), \\code{Y} (oracle target",
    "#'   phenotype), \\code{Y_obs} (measured target phenotype, \\code{NA} where",
    "#'   unmeasured), \\code{S} (synthetic phenotype) and \\code{GRM} (sparse",
    "#'   genetic relatedness matrix).",
    "#'",
    "#' @return A named list of precomputed step-1 quantities, to be passed as",
    "#'   \\code{step1_pars} to \\code{\\link{score_test_ObsG_hatY_multiply}}."),

  unObsG_hatY_typeIerror_step1 = c(
    "#' Null model for the synthetic phenotype among unmeasured individuals",
    "#'",
    "#' As \\code{\\link{ObsG_hatY_typeIerror_step1}}, but restricted to the",
    "#' individuals with no measured target phenotype. Fitting the two subsets",
    "#' separately makes it possible to check that the synthetic phenotype behaves",
    "#' the same way in both.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all}, \\code{Y}, \\code{Y_obs},",
    "#'   \\code{S} and \\code{GRM}; see \\code{\\link{ObsG_hatY_typeIerror_step1}}.",
    "#'",
    "#' @return A named list of precomputed step-1 quantities, to be passed as",
    "#'   \\code{step1_pars} to \\code{\\link{score_test_unObsG_hatY_multiply}}."),

  score_test_ObsG_hatY_multiply = c(
    "#' Score test of the synthetic phenotype among measured individuals",
    "#'",
    "#' Mixed-model score test for association between each variant and the",
    "#' synthetic phenotype, restricted to individuals who also have a measured",
    "#' target phenotype.",
    "#'",
    "#' @param g_matrix Numeric genotype matrix, individuals in rows and variants",
    "#'   in columns, over the full cohort.",
    "#' @param step1_pars The list returned by",
    "#'   \\code{\\link{ObsG_hatY_typeIerror_step1}}.",
    "#'",
    "#' @return A named list with the score statistic, \\eqn{-\\log_{10}} p-value,",
    "#'   effect estimate and its variance, one entry per variant."),

  score_test_unObsG_hatY_multiply = c(
    "#' Score test of the synthetic phenotype among unmeasured individuals",
    "#'",
    "#' Mixed-model score test for association between each variant and the",
    "#' synthetic phenotype, restricted to individuals with no measured target",
    "#' phenotype.",
    "#'",
    "#' @param g_matrix Numeric genotype matrix, individuals in rows and variants",
    "#'   in columns, over the full cohort.",
    "#' @param step1_pars The list returned by",
    "#'   \\code{\\link{unObsG_hatY_typeIerror_step1}}.",
    "#'",
    "#' @return A named list with the score statistic, \\eqn{-\\log_{10}} p-value,",
    "#'   effect estimate and its variance, one entry per variant."),

  ObsG_typeIerror_step1_lm = c(
    "#' Ordinary least squares null model for the measured target phenotype",
    "#'",
    "#' Fits the null model by ordinary least squares, ignoring relatedness, to the",
    "#' individuals with a measured target phenotype. Provides the reference",
    "#' against which the mixed-model analyses are compared: any inflation this",
    "#' shows but the mixed model does not is attributable to relatedness.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all}, \\code{Y}, \\code{Y_obs}",
    "#'   and \\code{S}. The genetic relatedness matrix is not used.",
    "#'",
    "#' @return A named list holding the residuals, residual variance estimate and",
    "#'   \\eqn{(X'X)^{-1}}, to be passed as \\code{step1_pars} to",
    "#'   \\code{\\link{score_test_ObsG_lm_multiply}}."),

  unObsG_hatY_typeIerror_step1_lm = c(
    "#' Ordinary least squares null model for the synthetic phenotype",
    "#'",
    "#' As \\code{\\link{ObsG_typeIerror_step1_lm}}, but fitted to the synthetic",
    "#' phenotype among the individuals with no measured target phenotype.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all}, \\code{Y}, \\code{Y_obs}",
    "#'   and \\code{S}. The genetic relatedness matrix is not used.",
    "#'",
    "#' @return A named list holding the residuals, residual variance estimate and",
    "#'   \\eqn{(X'X)^{-1}}, to be passed as \\code{step1_pars} to",
    "#'   \\code{\\link{score_test_unObsG_hatY_lm_multiply}}."),

  score_test_ObsG_lm_multiply = c(
    "#' Ordinary least squares score test, measured individuals",
    "#'",
    "#' Score test under ordinary linear regression, ignoring relatedness, for the",
    "#' measured target phenotype.",
    "#'",
    "#' @param g_matrix Numeric genotype matrix, individuals in rows and variants",
    "#'   in columns, over the full cohort.",
    "#' @param step1_pars The list returned by",
    "#'   \\code{\\link{ObsG_typeIerror_step1_lm}}.",
    "#'",
    "#' @return A named list with the score statistic, \\eqn{-\\log_{10}} p-value,",
    "#'   effect estimate and its variance, one entry per variant."),

  score_test_unObsG_hatY_lm_multiply = c(
    "#' Ordinary least squares score test, synthetic phenotype",
    "#'",
    "#' Score test under ordinary linear regression, ignoring relatedness, for the",
    "#' synthetic phenotype among individuals with no measured target phenotype.",
    "#'",
    "#' @param g_matrix Numeric genotype matrix, individuals in rows and variants",
    "#'   in columns, over the full cohort.",
    "#' @param step1_pars The list returned by",
    "#'   \\code{\\link{unObsG_hatY_typeIerror_step1_lm}}.",
    "#'",
    "#' @return A named list with the score statistic, \\eqn{-\\log_{10}} p-value,",
    "#'   effect estimate and its variance, one entry per variant."),

  SynSurrG_eachG_step1 = c(
    "#' SynPALM null model with the tested variant retained in the surrogate model",
    "#'",
    "#' As \\code{\\link{SynSurrG_ablation_estimate}}, except that the variant under",
    "#' test is kept as a covariate when the synthetic phenotype is regressed on the",
    "#' covariates. This is the appropriate step 1 when the prediction model that",
    "#' generated the synthetic phenotype already contains that variant, so that the",
    "#' variance components are not contaminated by the signal being tested. Because",
    "#' step 1 then depends on the variant, it must be refitted for each one, which",
    "#' is far more costly than the shared-step-1 route.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all}, \\code{Y}, \\code{Y_obs},",
    "#'   \\code{S} and \\code{GRM}.",
    "#' @param G Numeric genotype vector for the single variant to be retained as a",
    "#'   covariate in the surrogate model.",
    "#'",
    "#' @return A named list of precomputed step-1 quantities, to be passed as",
    "#'   \\code{step1_pars} to \\code{\\link{score_test_SynSurrG_single}}."),

  Oracle_ablation_estimate = c(
    "#' Oracle null model on independent samples",
    "#'",
    "#' Step 1 for the oracle analysis, which uses the complete target phenotype",
    "#' that would be available if nobody were missing a measurement. Restricted to",
    "#' one individual per relatedness block, so the covariance is diagonal and no",
    "#' genetic relatedness matrix is needed. Provides the upper bound on power in",
    "#' the ablation comparisons.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all} (covariates, columns prefixed",
    "#'   \"X.\"), \\code{Y} (complete target phenotype), \\code{Y_obs} and \\code{S}.",
    "#' @param independent_indices Integer vector of row indices giving one",
    "#'   individual per relatedness block, as returned by sampling within the",
    "#'   output of \\code{\\link{find_blocks_vectorized}}.",
    "#'",
    "#' @return A named list of precomputed step-1 quantities, to be passed as",
    "#'   \\code{step1_pars} to \\code{\\link{score_test_Oracle_multiply}}."),

  Obs_ablation_estimate = c(
    "#' Observed-only null model on independent samples",
    "#'",
    "#' Step 1 for the observed-only analysis, which discards the synthetic",
    "#' phenotype and uses just the individuals with a measured target phenotype.",
    "#' Restricted to one individual per relatedness block, so the covariance is",
    "#' diagonal and no genetic relatedness matrix is needed. Provides the baseline",
    "#' against which the power gain from the synthetic phenotype is measured.",
    "#'",
    "#' @param mydf A list with elements \\code{X_all} (covariates, columns prefixed",
    "#'   \"X.\"), \\code{Y}, \\code{Y_obs} (\\code{NA} where unmeasured) and \\code{S}.",
    "#' @param independent_indices Integer vector of row indices giving one",
    "#'   individual per relatedness block, as returned by sampling within the",
    "#'   output of \\code{\\link{find_blocks_vectorized}}.",
    "#'",
    "#' @return A named list of precomputed step-1 quantities, to be passed as",
    "#'   \\code{step1_pars} to \\code{\\link{score_test_Obs_multiply}}.")
)

## ---------------------------------------------------------------------------

lines <- readLines(TARGET, warn = FALSE)

def_line <- function(nm, lines) {
  hit <- grep(paste0("^", nm, "[ \t]*(<-|=)[ \t]*function[ \t]*\\("), lines)
  if (length(hit) != 1L) return(NA_integer_)
  hit
}

replaced <- character(0); skipped <- character(0)

## work from the bottom up so earlier line numbers stay valid
order_by_pos <- names(docs)[order(vapply(names(docs),
                                         function(n) def_line(n, lines),
                                         numeric(1)),
                                  decreasing = TRUE)]

for (nm in order_by_pos) {
  d <- def_line(nm, lines)
  if (is.na(d)) { skipped <- c(skipped, nm); next }

  ## extent of the existing roxygen block above the definition
  top <- d
  k <- d - 1L
  while (k >= 1L && grepl("^\\s*#'", lines[k])) { top <- k; k <- k - 1L }

  new_block <- c(docs[[nm]], "#' @export")

  if (top == d) {
    lines <- append(lines, new_block, after = d - 1L)   # no block present
  } else {
    lines <- c(lines[seq_len(top - 1L)], new_block, lines[d:length(lines)])
  }
  replaced <- c(replaced, nm)
}

writeLines(lines, TARGET)

cat("Documentation written for ", length(replaced), " functions.\n", sep = "")
if (length(skipped)) {
  cat("NOT FOUND in the file (check the names):\n")
  cat(paste0("  ", skipped, collapse = "\n"), "\n")
}

ok <- tryCatch({ parse(TARGET); TRUE },
               error = function(e) { cat("PARSE ERROR: ", conditionMessage(e), "\n"); FALSE })
cat("File still parses cleanly: ", ok, "\n", sep = "")

left <- grep("TODO: one line describing", readLines(TARGET, warn = FALSE))
cat("Remaining placeholder stubs: ", length(left), "\n", sep = "")

cat("\nNext: devtools::document() then devtools::check()\n")

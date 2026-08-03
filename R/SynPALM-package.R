#' SynPALM: Synthetic Phenotype Assisted Linear Mixed Models
#'
#' Robust and computationally scalable inference for proteome-wide GWAS when the
#' target phenotype is only partially observed. SynPALM jointly analyses observed
#' measurements and complete synthetic (machine-learning predicted) phenotypes
#' while accounting for cryptic relatedness and population structure via linear
#' mixed models.
#'
#' @import Matrix
#' @importFrom methods as
#' @importFrom parallel mclapply
#' @importFrom dplyr %>%
#' @importFrom stats complete.cases cov lm na.omit pchisq qnorm rbinom
#'   reformulate rnorm setNames var
#' @importFrom utils write.table globalVariables
#'
#' @keywords internal
"_PACKAGE"

## ---------------------------------------------------------------------------
## Many functions receive a list of precomputed quantities and expand it into
## the local frame with
##     list2env(step1_pars, envir = environment())
## The static code analyser in R CMD check cannot see bindings created that way
## and reports them as undefined globals. Declaring them here silences those
## notes. This list is the union of every name passed through list2env().
## ---------------------------------------------------------------------------
utils::globalVariables(c(
  "A22", "Atb", "Att", "B1", "B2", "B2_1", "B_mat", "Btt1",
  "L_cond", "L_cond_oracle",
  "SX", "Sigma11", "Sigma11_oracle", "Sigma12", "Sigma22",
  "Sigma12_Sigma22inv", "Sigma12_oracle_Sigma22inv",
  "V11", "X", "X_obs", "XtSX_inv", "XtX_inv", "Y",
  "bt2", "chol_Sigma22", "hatY", "hat_sigma2",
  "invS11_res", "inv_Sigma11", "inv_Sigma22",
  "n_obs", "n_unobs", "obs_protein_index", "residual",
  "unobs_protein_index"
))

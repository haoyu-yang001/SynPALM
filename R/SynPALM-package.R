#' SynPALM: Synthetic Phenotype Assisted Linear Mixed Models
#'
#' Robust and computationally scalable inference for proteome-wide GWAS when the
#' target phenotype is only partially observed. SynPALM jointly analyses observed
#' measurements and complete synthetic (machine-learning predicted) phenotypes
#' while accounting for cryptic relatedness and population structure via linear
#' mixed models.
#'
#' @section Pipeline:
#' \code{\link{synpalm_gwas}} runs a whole analysis of one protein:
#' \code{\link{synpalm_folds}} (relatedness-aware folds),
#' \code{\link{synpalm_predict}} (cross-fitted random-forest synthetic
#' phenotype), \code{\link{synpalm_null}} (null model) and
#' \code{\link{synpalm_scan}} (score tests). The steps can also be called one
#' at a time, e.g. to fit the null model once and scan chromosomes in
#' parallel jobs.
#'
#' @section Numerical stability:
#' The Haseman-Elston moment estimators of the variance components are
#' unconstrained and can be negative or imply an indefinite joint covariance.
#' The step-1 and ablation functions pass their estimates through
#' \code{\link{constrain_vc_bivariate}} or \code{\link{constrain_vc_univariate}}.
#' By default (\code{"none"}) the raw estimates are kept, with only the
#' residual variance of the protein floored at 0.01, as in the original
#' analysis; \code{"project"} instead projects them onto the set of valid
#' covariances. Either way \code{vc} records whether the estimates were valid.
#' Block-wise inverses (\code{\link{matrix_inv_block}},
#' \code{\link{matrix_inv_Amatrix}}) catch failures per block and regularise
#' blocks that are not positive definite.
#'
#' Every \code{*_step1} and \code{*_estimate} function also returns a
#' \code{vc} element holding the raw and the used variance components and the
#' projection flags, so a regularised protein can be identified afterwards.
#'
#' Behaviour is controlled with \code{options()}:
#' \describe{
#'   \item{\code{synsurrg.on_nonpsd}}{\code{"none"} (package default) keeps
#'     the raw estimates; \code{"project"} regularises invalid estimates;
#'     \code{"stop"} raises an error. \code{\link{synpalm_null}} sets it per
#'     call through its \code{vc_constraint} argument.}
#'   \item{\code{synsurrg.verbose}}{\code{TRUE} (default) prints variance-component
#'     and block-inversion diagnostics.}
#'   \item{\code{synsurrg.ridge_e}}{Residual-variance floor as a fraction of
#'     the mean squared residual. Default 0.01.}
#'   \item{\code{synsurrg.ncores}}{Cores for block-wise inversion. Default
#'     \code{NULL}, meaning \code{SLURM_CPUS_PER_TASK} if set, else 2.}
#' }
#'
#' @import Matrix
#' @importFrom methods as is
#' @importFrom parallel mclapply
#' @importFrom dplyr %>%
#' @importFrom stats complete.cases cov lm na.exclude na.omit pchisq qnorm rbinom
#'   reformulate residuals rnorm setNames var
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

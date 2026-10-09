## Package default for the variance-component constraint: "none", i.e. the raw
## Haseman-Elston estimates with only sigma_T2 floored at 0.01, as in the code
## used before 2026-08-22. A value the user has already set is left alone.
## See ?SynPALM-package, section "Numerical stability".
.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("synsurrg.on_nonpsd")))
    options(synsurrg.on_nonpsd = "none")
  invisible()
}

## ---------------------------------------------------------------------------
## Collects every number the Nature Research software submission checklist
## requires, and prints them in a form you can paste straight into README.md.
##
## Run from an R session that already has SynPALM installed:
##   source("tools/collect_checklist_numbers.R")
##
## Set MEASURE_CLEAN_INSTALL <- TRUE to also time a from-scratch install into a
## temporary library. That downloads and compiles every dependency and can take
## 10-30 minutes, so it is off by default. Run it once, on a normal
## desktop/laptop -- NOT on a cluster login node, since the checklist asks for
## timings on a "normal desktop computer".
## ---------------------------------------------------------------------------

MEASURE_CLEAN_INSTALL <- FALSE

deps <- c("Matrix", "stats", "Rcpp", "lattice", "data.table",
          "dplyr", "GMMAT", "ranger", "withr")

cat("\n==========================================================\n")
cat("1. SYSTEM REQUIREMENTS\n")
cat("==========================================================\n\n")

cat("R version:  ", R.version.string, "\n", sep = "")
si <- Sys.info()
cat("OS:         ", si[["sysname"]], " ", si[["release"]], "\n", sep = "")
cat("Machine:    ", si[["machine"]], "\n\n", sep = "")

vers <- vapply(deps, function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)),
                error = function(e) "NOT INSTALLED")
  v
}, character(1))

cat("--- paste this table into README.md ---\n\n")
cat("| Package | Version tested |\n|---|---|\n")
for (p in deps) cat("| ", p, " | ", vers[[p]], " |\n", sep = "")

cat("\n--- paste these lines into DESCRIPTION (version floors) ---\n\n")
cat("Depends:\n    R (>= ", paste(R.version$major, R.version$minor, sep = "."),
    ")\n", sep = "")
cat("Imports:\n")
imp <- setdiff(deps, "stats")
for (i in seq_along(imp)) {
  p <- imp[i]
  comma <- if (i < length(imp)) "," else ""
  cat("    ", p, " (>= ", vers[[p]], ")", comma, "\n", sep = "")
}
cat("    stats\n")

cat("\n==========================================================\n")
cat("2. MEMORY FOOTPRINT OF THE DEMO\n")
cat("==========================================================\n\n")

gc(reset = TRUE, full = TRUE)
invisible(capture.output(
  source(system.file("demo", "quickstart.R", package = "SynPALM"))
))
peak_mb <- sum(gc()[, "max used"] * c(8, 8) / 1e6)
cat("Peak R memory during demo (approx. MB): ", round(peak_mb), "\n", sep = "")
cat("(round up generously for the README; this excludes OS overhead)\n")

cat("\n==========================================================\n")
cat("3. DEMO OUTPUT AND RUN TIME\n")
cat("==========================================================\n\n")
cat("Running the seeded demo again with full output.\n")
cat("Copy EVERYTHING between the markers into the README\n")
cat("'Expected output' block, verbatim.\n\n")
cat("-------- BEGIN EXPECTED OUTPUT --------\n")
source(system.file("demo", "quickstart.R", package = "SynPALM"))
cat("--------- END EXPECTED OUTPUT ---------\n")

cat("\n==========================================================\n")
cat("4. INSTALL TIME\n")
cat("==========================================================\n\n")

if (isTRUE(MEASURE_CLEAN_INSTALL)) {
  tmplib <- file.path(tempdir(), paste0("synpalm_lib_", as.integer(Sys.time())))
  dir.create(tmplib, recursive = TRUE)
  cat("Installing into a clean library at:\n  ", tmplib, "\n", sep = "")
  cat("This will take a while. Do not interrupt.\n\n")

  t_clean <- system.time(
    withr::with_libpaths(c(tmplib), action = "prefix", {
      remotes::install_github("haoyu-yang001/SynPALM",
                              dependencies = TRUE,
                              upgrade = "never",
                              quiet = FALSE)
    })
  )
  cat("\nClean-library install, elapsed (minutes): ",
      round(t_clean[["elapsed"]] / 60, 1), "\n", sep = "")
  unlink(tmplib, recursive = TRUE)
} else {
  cat("SKIPPED. Set MEASURE_CLEAN_INSTALL <- TRUE at the top and re-run.\n")
}

cat("\nNow timing a reinstall with dependencies already present.\n\n")
t_warm <- system.time(
  remotes::install_github("haoyu-yang001/SynPALM",
                          dependencies = FALSE,
                          upgrade = "never",
                          quiet = TRUE)
)
cat("Install with dependencies present, elapsed (minutes): ",
    round(t_warm[["elapsed"]] / 60, 1),
    "  (= ", round(t_warm[["elapsed"]]), " seconds)\n", sep = "")

cat("\n==========================================================\n")
cat("DONE. Save this whole console transcript.\n")
cat("==========================================================\n")

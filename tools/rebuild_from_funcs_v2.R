## ---------------------------------------------------------------------------
## Rebuild the package's R code from funcs_v2.R, the code that produced the
## manuscript results.
##
## Function boundaries are determined by R's own parser via srcref, not by
## counting braces, so nested functions and multi-line calls are handled
## exactly as R sees them.
##
## Bodies are copied verbatim. Documentation is carried over from the previous
## package files where a matching function exists, otherwise a stub is written.
##
## The only deliberate change to any body is replacing a non-ASCII character
## inside a message() string, which R CMD check rejects and which prints as
## garbage on Windows. Nothing affecting a computed value is touched.
##
## Safe to re-run.
##
## Usage, from the package root:
##   source("tools/rebuild_from_funcs_v2.R")
## ---------------------------------------------------------------------------

SRC <- path.expand("~/Desktop/funcs_v2.R")
OUT <- "R/synpalm_functions.R"

stopifnot("funcs_v2.R not found on the Desktop" = file.exists(SRC))

cand <- c("R/core_functions.R", "R/utils.R",
          "tools/old_R_backup/core_functions.R", "tools/old_R_backup/utils.R")
OLD_FILES <- cand[file.exists(cand)]
OLD_FILES <- OLD_FILES[!duplicated(basename(OLD_FILES))]
cat("Reusing documentation from:\n")
cat(paste0("  ", OLD_FILES, collapse = "\n"), "\n\n")

## --- helpers ---------------------------------------------------------------

## Is this top-level expression a `name <- function(...)` assignment?
is_fun_def <- function(e) {
  is.call(e) && length(e) >= 3 &&
    as.character(e[[1]])[1] %in% c("<-", "=", "<<-") &&
    is.symbol(e[[2]]) &&
    is.call(e[[3]]) && identical(as.character(e[[3]][[1]])[1], "function")
}

## Every top-level function definition in a file, with its exact source text
## and starting line, taken from R's parser.
parse_defs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  p <- parse(path, keep.source = TRUE)
  refs <- attr(p, "srcref")
  out <- list()
  for (i in seq_along(p)) {
    e <- p[[i]]
    if (!is_fun_def(e)) next
    nm <- as.character(e[[2]])
    r  <- refs[[i]]
    out[[nm]] <- list(
      name  = nm,
      start = r[1L],
      code  = lines[r[1L]:r[3L]],
      args  = { a <- names(e[[3]][[2]]); if (is.null(a)) character(0) else a[nzchar(a)] }
    )
  }
  list(defs = out, lines = lines)
}

roxygen_above <- function(lines, start) {
  k <- start - 1L; out <- character(0)
  while (k >= 1L && grepl("^\\s*#'", lines[k])) { out <- c(lines[k], out); k <- k - 1L }
  out
}

make_stub <- function(d) {
  c(paste0("#' ", d$name),
    "#'",
    paste0("#' TODO: one line describing what ", d$name, "() does."),
    "#'",
    if (length(d$args)) paste0("#' @param ", d$args, " TODO: describe.") else character(0),
    "#'",
    "#' @return TODO: describe the return value.",
    "#' @export")
}

ensure_export <- function(rox) if (any(grepl("@export", rox))) rox else c(rox, "#' @export")

## --- read both sides -------------------------------------------------------

new <- parse_defs(SRC)
new_defs <- new$defs

old_rox <- list()
for (f in OLD_FILES) {
  o <- parse_defs(f)
  for (nm in names(o$defs)) {
    rox <- roxygen_above(o$lines, o$defs[[nm]]$start)
    if (length(rox)) old_rox[[nm]] <- rox
  }
}

cat("Top-level functions in funcs_v2.R : ", length(new_defs), "\n", sep = "")
cat("Roxygen blocks available to reuse : ", length(old_rox), "\n\n", sep = "")

## --- assemble --------------------------------------------------------------

header <- c(
  "## ---------------------------------------------------------------------------",
  "## SynPALM: function definitions.",
  "##",
  "## The code below is taken verbatim from the analysis source used to produce",
  "## the manuscript results, so that the package and the published analysis are",
  "## guaranteed to agree.",
  "## ---------------------------------------------------------------------------",
  ""
)

body_out <- character(0); reused <- 0L; stubbed <- character(0)

for (nm in names(new_defs)) {
  d <- new_defs[[nm]]
  code <- gsub("\u2705 ", "Done: ", d$code, fixed = TRUE)
  code <- gsub("\u2705",  "Done:",  code,   fixed = TRUE)

  if (!is.null(old_rox[[nm]])) {
    rox <- ensure_export(old_rox[[nm]]); reused <- reused + 1L
  } else {
    rox <- make_stub(d); stubbed <- c(stubbed, nm)
  }
  body_out <- c(body_out, rox, code, "")
}

writeLines(c(header, body_out), OUT)

## --- verify ----------------------------------------------------------------

cat("Wrote ", OUT, "  (", length(c(header, body_out)), " lines)\n", sep = "")
cat("Roxygen reused : ", reused, "\n", sep = "")
cat("Stubs written  : ", length(stubbed), "\n\n", sep = "")

ok <- tryCatch({ parse(OUT); TRUE },
               error = function(e) { cat("PARSE ERROR: ", conditionMessage(e), "\n"); FALSE })
cat("File parses cleanly            : ", ok, "\n", sep = "")

if (ok) {
  chk <- parse(OUT, keep.source = TRUE)
  kinds <- vapply(as.list(chk), is_fun_def, logical(1))
  cat("Top-level expressions          : ", length(chk), "\n", sep = "")
  cat("Of which function definitions  : ", sum(kinds), "\n", sep = "")
  cat("Stray non-function top level   : ", sum(!kinds),
      if (sum(!kinds) == 0) "  (good)" else "  <-- PROBLEM", "\n", sep = "")
  if (any(!kinds)) {
    cat("\n  stray expressions:\n")
    for (i in which(!kinds))
      cat("    ", substr(paste(deparse(chk[[i]]), collapse = " "), 1, 90), "\n", sep = "")
  }
  written <- vapply(as.list(chk)[kinds], function(e) as.character(e[[2]]), character(1))
  cat("Duplicated names               : ", any(duplicated(written)), "\n", sep = "")
  miss <- setdiff(names(new_defs), written)
  cat("Missing vs funcs_v2.R          : ",
      if (length(miss)) paste(miss, collapse = ", ") else "none", "\n", sep = "")
}

if (length(stubbed)) {
  cat("\nFunctions with placeholder docs (need a one-line description):\n")
  cat(paste0("  ", stubbed, collapse = "\n"), "\n")
}

## --- retire the old files --------------------------------------------------

bak <- file.path("tools", "old_R_backup")
dir.create(bak, recursive = TRUE, showWarnings = FALSE)
for (f in c("R/core_functions.R", "R/utils.R")) {
  if (file.exists(f)) {
    file.copy(f, file.path(bak, basename(f)), overwrite = TRUE)
    file.remove(f)
    cat("moved ", f, " -> ", file.path(bak, basename(f)), "\n", sep = "")
  }
}

cat("\nR/ now contains:\n"); print(list.files("R"))
cat("\nNext: devtools::document() then devtools::check()\n")

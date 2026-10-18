## Type I error tables and QQ plots from typeIerror_rev1.R (Table 1, Fig 2B,
## Supp Figs S2-S4). Streams the array files: exact exceedance counts at each
## alpha, binned -log10 p for the bulk of the QQ plot and the exact values of
## the tail (-log10 p > 4).
##   Rscript analyze_t1e_rev1.R
## Writes simulation_code/results_rev1/figs/: table_typeIerror_rev1.csv (all r, pi,
## rho, methods), table1_rev1.csv (sib pairs, Syn-PALM and POP-GWAS at 1e-3,
## 1e-5, 1e-7 as in Table 1) and qq_GRMr_<r>_SynG_rev1.{pdf,png}.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
D   <- Sys.getenv("T1E_DIR", "simulation_code/results_rev1/typeIerror/")   # typeIerror_split/ for the split-pair design
FIG <- Sys.getenv("T1E_FIG", "simulation_code/results_rev1/figs/"); dir.create(FIG, showWarnings = FALSE)
ALPHA <- 10^-(2:7)
BRK   <- seq(0, 4, by = 0.002)
files <- list.files(D, "_array_\\d+\\.rds$", full.names = TRUE)
meta  <- data.table(f = files,
  GRMr = as.numeric(sub("^GRMr_([0-9.]+)_.*", "\\1", basename(files))),
  miss = as.numeric(sub(".*_missing_([0-9.]+)_array.*", "\\1", basename(files))),
  array = as.integer(sub(".*_array_(\\d+)\\.rds", "\\1", basename(files))))
meta <- meta[array <= as.integer(Sys.getenv("T1E_MAX_ARRAY", "500"))]

## T1E_FROM_RDS=1: reuse the counts saved by a previous run (t1e_accumulated_rev1.rds)
## instead of streaming the array files again (tables and QQ plots only; same output)
from_rds <- Sys.getenv("T1E_FROM_RDS") == "1"
acc <- if (from_rds) readRDS(paste0(FIG, "t1e_accumulated_rev1.rds")) else list()   # key -> list(n, exceed counts, bin counts, tail values)
for (i in seq_len(nrow(meta) * !from_rds)) {
  x <- readRDS(meta$f[i])
  for (rn in names(x)) {
    m <- x[[rn]]
    nl <- list(SynSurrG = m[, "SynSurrG_negative_log10_pvalues"])
    if ("ObsG_negative_log10_pvalues" %in% colnames(m)) nl$ObsG <- m[, "ObsG_negative_log10_pvalues"]
    if ("POPGWAS_p" %in% colnames(m)) nl$POPGWAS <- -log10(m[, "POPGWAS_p"])
    for (me in names(nl)) {
      v <- nl[[me]]; v <- v[is.finite(v)]
      key <- paste(meta$GRMr[i], meta$miss[i], sub("rho_", "", rn), me, sep = "|")
      a <- acc[[key]]
      if (is.null(a)) a <- list(n = 0, ex = numeric(length(ALPHA)), bins = numeric(length(BRK) - 1), tail = numeric(0), n_files = 0)
      a$n <- a$n + length(v)
      a$ex <- a$ex + sapply(ALPHA, function(al) sum(v >= -log10(al)))
      a$bins <- a$bins + tabulate(findInterval(v[v <= 4], BRK, rightmost.closed = TRUE), length(BRK) - 1)
      a$tail <- c(a$tail, v[v > 4]); a$n_files <- a$n_files + 1
      acc[[key]] <- a
    }
  }
}
if (!from_rds) saveRDS(acc, paste0(FIG, "t1e_accumulated_rev1.rds"))   # for re-plotting without re-reading the arrays
k <- tstrsplit(names(acc), "|", fixed = TRUE)
tab <- data.table(GRMr = as.numeric(k[[1]]), miss = as.numeric(k[[2]]), rho = as.numeric(k[[3]]), method = k[[4]],
                  n_arrays = sapply(acc, `[[`, "n_files"), n_snps = sapply(acc, `[[`, "n"),
                  lambda_gc = NA_real_)
for (j in seq_along(ALPHA)) tab[, (formatC(ALPHA[j], format = "e", digits = 0)) := sapply(acc, function(a) a$ex[j] / a$n)]
## lambda_GC from the binned distribution (median of chi2_1)
## (linear interpolation within the bin that contains the median, so that lambda is
## not limited to the resolution of the 0.002 bins)
tab$lambda_gc <- sapply(acc, function(a) {
  above <- length(a$tail) + c(0, cumsum(rev(a$bins)))      # SNPs with -log10 p above each upper edge
  up <- rev(BRK)                                              # upper edges, from 4 down to 0
  k <- which(above >= a$n / 2)[1]                             # median lies in bin (up[k], up[k-1])
  f <- (a$n / 2 - above[k - 1]) / (above[k] - above[k - 1])
  med_nlp <- up[k - 1] - f * (up[k - 1] - up[k])
  qchisq(10^-med_nlp, 1, lower.tail = FALSE) / qchisq(0.5, 1) })
setorder(tab, GRMr, method, rho, miss)
fwrite(tab, paste0(FIG, "table_typeIerror_rev1.csv"))
t1 <- dcast(tab[GRMr == 0.5 & method %in% c("SynSurrG", "POPGWAS")], miss + rho ~ method,
            value.var = c("1e-03", "1e-05", "1e-07"))
fwrite(t1, paste0(FIG, "table1_rev1.csv"))
print(tab[, .(GRMr, method, miss, rho, n_arrays, n_snps, lambda_gc, `1e-03`, `1e-05`, `1e-07`)], digits = 3)

## QQ plots (Syn-PALM), one figure per r, facets by missing rate, colour by rho
qq_points <- function(a) {
  N <- a$n; tl <- sort(a$tail, decreasing = TRUE); nt <- length(tl)
  e_t <- -log10((seq_len(nt) - 0.5) / N)
  ## bulk: 3000 ranks spread over the rest, observed value from the bins
  cs <- cumsum(rev(a$bins)); mid <- rev((BRK[-1] + BRK[-length(BRK)]) / 2)
  rk <- unique(round(exp(seq(log(nt + 1), log(N), length.out = 3000))))
  o_b <- mid[findInterval(rk - nt, c(0, cs), left.open = TRUE)]
  rr <- c(seq_len(nt), rk)
  data.table(expected = -log10((rr - 0.5) / N), observed = c(tl, o_b),
             clower = -log10(qbeta(0.025, rr, N - rr + 1)), cupper = -log10(qbeta(0.975, rr, N - rr + 1)))
}
## common look: theme_rev1.R (rho light to dark blue, red identity line, grey 95% band
## with dashed borders from the rho = 0 ranks)
source("simulation_code/rev1/theme_rev1.R")
for (g in sort(unique(tab$GRMr))) {
  ks <- names(acc)[startsWith(names(acc), paste0(g, "|")) & endsWith(names(acc), "|SynSurrG")]
  if (!length(ks)) next
  pts <- rbindlist(lapply(ks, function(kk) { z <- strsplit(kk, "|", fixed = TRUE)[[1]]
    cbind(miss = as.numeric(z[2]), rho = as.numeric(z[3]), qq_points(acc[[kk]])) }))
  pts[, rho_f := factor(rho, sort(unique(rho)))]
  pts[, panel := factor(miss_lab(miss), miss_lab(sort(unique(miss))))]
  p <- ggplot(pts, aes(expected, observed)) +
    qq_ref_layers(pts[rho == 0]) +
    geom_point(aes(colour = rho_f), size = 1.1, stroke = 0) +
    facet_wrap(~ panel, nrow = 2) +
    scale_colour_manual(name = expression(rho), values = rho_cols,
                        guide = guide_legend(override.aes = list(size = 2.6))) +
    scale_x_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.03)) +
    labs(x = qq_xlab, y = qq_ylab) +
    theme_rev1(12, legend = if (uniqueN(pts$panel) == 5) "inside" else "right") +   # inside: the empty 6th cell
    theme(legend.position.inside = c(0.835, 0.24), legend.justification = c(0.5, 0.5))
  ggsave(sprintf("%sqq_GRMr_%s_SynG_rev1.pdf", FIG, g), p, height = 5, width = 7.2, units = "in", device = cairo_pdf)
  ggsave(sprintf("%sqq_GRMr_%s_SynG_rev1.png", FIG, g), p, height = 5, width = 7.2, units = "in", dpi = 300)
}
cat("wrote tables and QQ plots to", FIG, "\n")

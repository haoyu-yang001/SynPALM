## Table 1 (revision, 2026-10-19): ratio of the empirical type I error to the nominal
## level alpha (1 = exact calibration) of Syn-PALM and POP-GWAS for sib pairs, at
## alpha = 1e-5 and 1e-7, for the two labelling designs in one table:
##   whole: each pair labelled or unlabelled as a whole (typeIerror/, 1e8 null SNPs)
##   split: individuals labelled at random, relatives split (typeIerror_split/, 1e7 null
##          SNPs; 1e-7 is not estimable, about one expected exceedance)
## Reads the exceedance counts saved by analyze_t1e_rev1.R (t1e_accumulated_rev1.rds).
## Bold: above 1 at one-sided exact Poisson P < 0.01. Rates and lambda_GC of the split
## design stay in Supplementary Table S6 (tab:sim_split).
##   Rscript table1_merged_rev1.R
## Writes simulation_code/results_rev1/figs/table1_merged_rev1.{csv,tex}
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
ALPHA <- 10^-(2:7)
rd <- function(p, des) {
  acc <- readRDS(p); k <- tstrsplit(names(acc), "|", fixed = TRUE)
  d <- data.table(des = des, GRMr = as.numeric(k[[1]]), miss = as.numeric(k[[2]]), rho = as.numeric(k[[3]]),
                  method = k[[4]], n = sapply(acc, `[[`, "n"))
  for (a in c(1e-5, 1e-7)) {
    j <- which(abs(log10(ALPHA) - log10(a)) < 1e-9)
    cnt <- sapply(acc, function(z) z$ex[j])
    d[, (sprintf("ratio_%g", a)) := cnt / (n * a)]
    d[, (sprintf("p_%g", a)) := ppois(cnt - 1, n * a, lower.tail = FALSE)]   # P(X >= cnt)
    d[, (sprintf("count_%g", a)) := cnt]
  }
  d[GRMr == 0.5 & method %in% c("SynSurrG", "POPGWAS")]
}
x <- rbind(rd("simulation_code/results_rev1/figs/t1e_accumulated_rev1.rds", "whole"),
           rd("simulation_code/results_rev1/figs_split/t1e_accumulated_rev1.rds", "split"))
fwrite(x, "simulation_code/results_rev1/figs/table1_merged_rev1.csv")
f <- function(dd, me, m, r, a) {
  z <- x[des == dd & method == me & miss == m & rho == r]
  if (!nrow(z) || (dd == "split" && a == 1e-7)) return(NULL)
  v <- sprintf(if (a == 1e-7) "%.1f" else "%.2f", z[[sprintf("ratio_%g", a)]])   # about 10 expected exceedances at 1e-7
  if (z[[sprintf("p_%g", a)]] < 0.01) v <- sprintf("\\textbf{%s}", v)
  v
}
cell <- function(...) { v <- f(...); if (is.null(v)) "--" else v }
rows <- character(0)
for (r in c(0, 0.25, 0.5, 0.75)) {
  for (m in c(0, 0.5, 0.75)) {
    rows <- c(rows, paste(m, r, cell("whole", "SynSurrG", m, r, 1e-5), cell("whole", "SynSurrG", m, r, 1e-7), "",
                          cell("whole", "POPGWAS", m, r, 1e-5), cell("whole", "POPGWAS", m, r, 1e-7), "",
                          cell("split", "SynSurrG", m, r, 1e-5), cell("split", "POPGWAS", m, r, 1e-5), sep = " & "))
    rows[length(rows)] <- paste0(rows[length(rows)], if (m == 0.75 && r < 0.75) " \\\\ \\hline" else " \\\\")
  }
}
writeLines(rows, "simulation_code/results_rev1/figs/table1_merged_rev1.tex")
cat(rows, sep = "\n")
print(x[, .(des, method, miss, rho, n, ratio_1e.05 = round(`ratio_1e-05`, 3), p5 = signif(`p_1e-05`, 2),
            ratio_1e.07 = round(`ratio_1e-07`, 2), c7 = `count_1e-07`)][order(des, method, rho, miss)])

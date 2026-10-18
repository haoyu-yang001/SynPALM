## Fig 2B (revision): QQ plots of Syn-PALM for sib pairs at missing rates 0, 50% and
## 90% in one row (layout of the submitted figure), from the accumulated counts that
## analyze_t1e_rev1.R saves (t1e_accumulated_rev1.rds), with its binned QQ points.
##   Rscript simulation_code/rev1/plot_fig2b_rev1.R
## Writes simulation_code/results_rev1/figs/fig2b_qq_siblings_rev1.{pdf,png}
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
FIG <- "simulation_code/results_rev1/figs/"
acc <- readRDS(paste0(FIG, "t1e_accumulated_rev1.rds"))
BRK <- seq(0, 4, by = 0.002)
qq_points <- function(a) {          # as in analyze_t1e_rev1.R
  N <- a$n; tl <- sort(a$tail, decreasing = TRUE); nt <- length(tl)
  cs <- cumsum(rev(a$bins)); mid <- rev((BRK[-1] + BRK[-length(BRK)]) / 2)
  rk <- unique(round(exp(seq(log(nt + 1), log(N), length.out = 3000))))
  o_b <- mid[findInterval(rk - nt, c(0, cs), left.open = TRUE)]
  rr <- c(seq_len(nt), rk)
  data.table(expected = -log10((rr - 0.5) / N), observed = c(tl, o_b),
             clower = -log10(qbeta(0.025, rr, N - rr + 1)), cupper = -log10(qbeta(0.975, rr, N - rr + 1)))
}
ks <- names(acc)[startsWith(names(acc), "0.5|") & endsWith(names(acc), "|SynSurrG")]
pts <- rbindlist(lapply(ks, function(kk) { z <- strsplit(kk, "|", fixed = TRUE)[[1]]
  cbind(miss = as.numeric(z[2]), rho = as.numeric(z[3]), qq_points(acc[[kk]])) }))
pts <- pts[miss %in% c(0, 0.5, 0.9)]
## common look: theme_rev1.R (rho light to dark blue, red identity line, grey 95% band
## with dashed borders from the rho = 0 ranks, as before)
source("simulation_code/rev1/theme_rev1.R")
pts[, rho_f := factor(rho, sort(unique(rho)))]
pts[, panel := factor(miss_lab(miss), miss_lab(sort(unique(miss))))]
p <- ggplot(pts, aes(expected, observed)) +
  qq_ref_layers(pts[rho == 0]) +
  geom_point(aes(colour = rho_f), size = 1.1, stroke = 0) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(name = expression(rho), values = rho_cols,
                      guide = guide_legend(override.aes = list(size = 2.4))) +
  scale_x_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = scales::breaks_width(2), expand = expansion(mult = 0.03)) +
  labs(x = qq_xlab, y = qq_ylab) +
  theme_rev1(12, legend = "right")
ggsave(sprintf("%sfig2b_qq_siblings_rev1.pdf", FIG), p, height = 2.45, width = 7, units = "in", device = cairo_pdf)
ggsave(sprintf("%sfig2b_qq_siblings_rev1.png", FIG), p, height = 2.45, width = 7, units = "in", dpi = 300)

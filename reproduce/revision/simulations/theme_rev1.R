## Common look of the simulation figures of the revision (sourced by the plotting scripts):
## one theme, one method palette and one rho palette, so that all panels match.
suppressPackageStartupMessages(library(ggplot2))

theme_rev1 <- function(base_size = 11, legend = "bottom") {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.4),
          strip.text = element_text(face = "bold", size = rel(0.95), margin = margin(3, 3, 3, 3)),
          axis.title = element_text(size = rel(1)),
          axis.text = element_text(colour = "grey15", size = rel(0.85)),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
          plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.05)),
          plot.subtitle = element_text(hjust = 0.5, size = rel(0.9), colour = "grey25"),
          legend.position = legend,
          legend.title = element_text(size = rel(0.9), face = "bold"),
          legend.text = element_text(size = rel(0.85)),
          legend.key = element_blank(), legend.background = element_blank(),
          legend.key.width = unit(1.4, "lines"),
          legend.margin = margin(0, 4, 0, 4), legend.box.spacing = unit(4, "pt"),
          plot.margin = margin(5, 8, 5, 5))
}

## methods: Syn-PALM blue, Obs-GWAS orange, SynSurr green, REML / BOLT-REML for the VC figure
method_cols <- c("Syn-PALM" = "#1F5FA8", "Obs-GWAS" = "#E07B1A", "SynSurr" = "#2E9E5B",
                 "MoM" = "#1F5FA8", "REML" = "#2E9E5B", "BOLT-REML" = "#E07B1A",
                 "POP-GWAS" = "#8E44AD", "Oracle-GWAS" = "grey35", "Ind-GWAS" = "#4D5A66")
## rho: light to dark blue
rho_cols <- c("0" = "#9ECAE1", "0.25" = "#5A9FD0", "0.5" = "#2166AC", "0.75" = "#08306B")
rho_lty  <- c("0" = "dotted", "0.25" = "dotdash", "0.5" = "dashed", "0.75" = "solid")
miss_lab <- function(m) sprintf("Missing rate %g%%", 100 * as.numeric(m))

## ---- helpers shared by the restyled plotting scripts (2026-10-16) -------------------
method_shapes <- c("Syn-PALM" = 16, "Obs-GWAS" = 17, "SynSurr" = 15, "Ind-GWAS" = 18,
                   "POP-GWAS" = 25, "Oracle-GWAS" = 4)
## missing rate: light to dark blue (Syn-PALM at increasing missingness)
miss_cols <- c("0" = "#9ECAE1", "0.25" = "#6BAED6", "0.5" = "#3182BD", "0.75" = "#1C5C9E", "0.9" = "#08306B")
## relatedness of the pairs: light to dark purple
kin_lab  <- c(cousin = "Cousins", halfsiblings = "Half-siblings", siblings = "Siblings")
kin_cols <- c("Cousins" = "#D1C4E9", "Half-siblings" = "#9575CD", "Siblings" = "#4527A0")
## facet labeller for rho (Greek letter as text, so that the strip is bold like the others):
## facet_wrap(~ rho, labeller = rho_labeller)
rho_labeller <- as_labeller(function(r) sprintf("\u03c1 = %g", as.numeric(r)))
## plotmath for a significance level, e.g. alpha_expr(5e-8) -> alpha == 5 %*% 10^-8
alpha_parts <- function(a) { e <- floor(log10(a) + 1e-9); list(m = signif(a / 10^e, 3), e = e) }
alpha_expr <- function(a) { z <- alpha_parts(a)
  if (z$m == 1) bquote(alpha == 10^.(z$e)) else bquote(alpha == .(z$m) %*% 10^.(z$e)) }
power_lab <- function(a) { z <- alpha_parts(a)
  if (z$m == 1) bquote("Power at" ~ alpha == 10^.(z$e)) else bquote("Power at" ~ alpha == .(z$m) %*% 10^.(z$e)) }
## axis titles of QQ plots
qq_xlab <- quote("Expected" ~ "\u2212" * log[10] * italic(P))
qq_ylab <- quote("Observed" ~ "\u2212" * log[10] * italic(P))
## reference layers of a QQ plot: 95% band (grey, dashed black borders) from `band`
## (columns expected, clower, cupper; NULL for none) and the red identity line
qq_ref_layers <- function(band = NULL) {
  l <- list()
  if (!is.null(band)) l <- list(
    geom_ribbon(data = band, aes(x = expected, ymin = clower, ymax = cupper), inherit.aes = FALSE,
                fill = "grey88", colour = NA),
    geom_line(data = band, aes(x = expected, y = cupper), inherit.aes = FALSE,
              linetype = "dashed", colour = "grey15", linewidth = 0.35),
    geom_line(data = band, aes(x = expected, y = clower), inherit.aes = FALSE,
              linetype = "dashed", colour = "grey15", linewidth = 0.35))
  c(l, list(geom_abline(intercept = 0, slope = 1, colour = "#D62728", linewidth = 0.5)))
}
## y breaks of a ratio axis that include 1
## (1 is not added when it would crowd a neighbouring break)
ratio_breaks_rev1 <- function(limits) {
  b <- scales::pretty_breaks(n = 5)(limits)
  if (min(abs(b - 1)) > 0.08 * diff(range(limits, na.rm = TRUE))) b <- c(b, 1)
  sort(unique(b))
}

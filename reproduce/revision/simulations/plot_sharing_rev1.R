## Figure for S3/S6 (genetic vs residual sharing): ratio of the mean model-based
## SE and power ratio at 5e-8, Syn-PALM vs Obs-GWAS, against rho_g, by pi and rho.
##   Rscript simulation_code/rev1/plot_sharing_rev1.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
x <- fread("simulation_code/results_rev1/summary_S36.csv")[method == "SynSurrG" & !is.na(rho_g)]
d <- rbind(x[, .(rho_tot, rho_g, miss, measure = "SE ratio (Syn-PALM / Obs-GWAS)", value = se_ratio_vs_ObsG)],
           x[, .(rho_tot, rho_g, miss, measure = "Power ratio at 5e-8", value = power_ratio_vs_ObsG)])
## common look: theme_rev1.R (missing rate light to dark blue); power ratio on top, SE ratio
## below (as in the caption), the row labels on the left as y titles
source("simulation_code/rev1/theme_rev1.R")
d[, `:=`(pi = factor(as.character(miss), names(miss_cols)), panel = sprintf("\u03c1 = %g", rho_tot))]
d[, measure := factor(measure, c("Power ratio at 5e-8", "SE ratio (Syn-PALM / Obs-GWAS)"),
                      c('"Power ratio"', '"SE ratio"'))]
p <- ggplot(d, aes(rho_g, value, colour = pi, shape = pi)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.45) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_grid(measure ~ panel, scales = "free_y", switch = "y",
             labeller = labeller(measure = label_parsed, panel = label_value)) +
  scale_colour_manual(name = "Missing rate", values = miss_cols, labels = function(m) sprintf("%g%%", 100 * as.numeric(m)), drop = TRUE) +
  scale_shape_manual(name = "Missing rate", values = c(`0` = 16, `0.25` = 18, `0.5` = 15, `0.75` = 18, `0.9` = 17),
                     labels = function(m) sprintf("%g%%", 100 * as.numeric(m)), drop = TRUE) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1), labels = function(x) format(x, drop0trailing = TRUE)) +
  labs(x = expression("Genetic correlation" ~ rho[g] ~ "(residual correlation adjusted to keep" ~ rho * ")"), y = NULL,
       subtitle = quote("Syn-PALM relative to Obs-GWAS; power at" ~ alpha == 5 %*% 10^-8)) +
  theme_rev1(12, legend = "top") +
  theme(strip.placement = "outside", strip.background.y = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = "plain", size = rel(1)),
        plot.subtitle = element_text(hjust = 0, size = rel(0.95)), panel.spacing.x = unit(10, "pt"))
dir.create("simulation_code/results_rev1/figs", showWarnings = FALSE)
ggsave("simulation_code/results_rev1/figs/S36_sharing.pdf", p, width = 6.6, height = 5, device = cairo_pdf)
ggsave("simulation_code/results_rev1/figs/S36_sharing.png", p, width = 6.6, height = 5, dpi = 300)

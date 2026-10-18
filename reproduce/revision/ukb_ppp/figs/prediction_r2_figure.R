## Out-of-fold prediction accuracy of the cross-fitted random-forest synthetic
## phenotypes across the 2,916 UKB-PPP proteins (no PRS / no PCs): R^2 = rho_oof^2
## (signed rho_oof < 0 counted as R^2 = 0 in the bands).
##   Rscript code/post_cf3/figs/prediction_r2_figure.R
## Writes results/figs_cf_noPRS_noPC/prediction_R2_across_proteins.{pdf,png} and
## results/table_results_cf_noPRS_noPC/prediction_R2_summary.txt
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(cowplot) })
setwd("<PROJECT_DIR>")
x <- as.data.table(readRDS("results/rho_oof_cf_noPRS_noPC.rds"))
x[, R2 := pmax(rho_oof, 0)^2]
x[, band := cut(R2, c(-Inf, 0.01, 0.1, 0.25, Inf), right = FALSE,
                labels = c("poor", "weak", "moderate", "good"))]
setorder(x, -R2); x[, rank := .I]
x[, label := toupper(gsub("_", "/", protein))]
cols <- c(poor = "#BDBDBD", weak = "#9ECAE1", moderate = "#4292C6", good = "#08306B")
labs_band <- c(poor = expression("Poor (" * R^2 * " < 0.01)"), weak = expression("Weak (" * R^2 * " = 0.01-0.1)"),
               moderate = expression("Moderate (" * R^2 * " = 0.1-0.25)"), good = expression("Good (" * R^2 * " > 0.25)"))
## examples: the five best, three around the median and three of the poorest with a known strong cis-pQTL
ex <- c(x[1:5, protein], x[round(.N / 2), protein],
        intersect(c("folr3", "lilrb5", "abo"), x$protein))
pa <- ggplot(x, aes(rank, R2, colour = band)) + geom_point(size = 0.6) +
  geom_hline(yintercept = median(x$R2), linetype = 2, colour = "grey40") +
  annotate("text", x = nrow(x) * 0.98, y = median(x$R2), label = sprintf("median~R^2 == %.3f", median(x$R2)), parse = TRUE,
           hjust = 1, vjust = -0.5, size = 3, colour = "grey30") +
  geom_text_repel(data = x[protein %in% ex], aes(label = label), colour = "black", size = 2.8, fontface = "italic",
                  min.segment.length = 0, box.padding = 0.4, nudge_x = 250, direction = "y", seed = 1, max.overlaps = Inf) +
  scale_colour_manual(values = cols, labels = labs_band, name = NULL) +
  labs(x = "Proteins ranked by prediction accuracy", y = expression("Out-of-fold prediction " * R^2)) +
  theme_classic() + theme(legend.position = c(0.72, 0.8), legend.background = element_blank())
pb <- ggplot(x, aes(R2, fill = band)) + geom_histogram(binwidth = 0.01, boundary = 0, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(x = expression("Out-of-fold prediction " * R^2), y = "Number of proteins") + theme_classic()
p <- plot_grid(pa, pb, labels = c("a", "b"), rel_widths = c(1.35, 1), label_size = 13)
for (ext in c("pdf", "png")) ggsave(sprintf("results/figs_cf_noPRS_noPC/prediction_R2_across_proteins.%s", ext), p, width = 10, height = 4.2, dpi = 300)
sink("results/table_results_cf_noPRS_noPC/prediction_R2_summary.txt", split = TRUE)
cat(sprintf("proteins %d; rho_oof median %.3f (Q1 %.3f, Q3 %.3f); R2 median %.3f, max %.3f (%s)\n", nrow(x), median(x$rho_oof),
            quantile(x$rho_oof, .25), quantile(x$rho_oof, .75), median(x$R2), max(x$R2), x$label[1]))
print(x[, .(proteins = .N, share = round(.N / nrow(x), 3)), by = band][order(band)])
cat("negative rho_oof:", x[rho_oof < 0, .N], "\n")
print(x[protein %in% ex, .(label, rho_oof = round(rho_oof, 3), R2 = round(R2, 3), rank)])
sink()

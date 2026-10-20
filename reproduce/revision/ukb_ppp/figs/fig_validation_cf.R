## Figure 6 (revision, 2026-10-19): validation of the UKB-PPP primary associations, the
## results of the "external validation" and "internal validation" paragraphs, in a compact
## 2 x 3 layout (the per-protein shares are quoted in the text).
##   A  share of the primary associations reported in the GWAS Catalog and with P < 5e-8
##      in the combined UKB-PPP analysis of Sun et al. (N = 51,507), pooled, for Obs-GWAS
##      and for the Syn-PALM associations shared with Obs-GWAS or additional
##   B  additional Syn-PALM associations: Obs-GWAS -log10 P against Syn-PALM -log10 P at
##      the Syn-PALM lead variant, coloured by P < 5e-8 in Sun et al.
##   C  additional associations: Obs-GWAS against Syn-PALM effect estimate
##   D  additional associations: ratio of the standard errors (Syn-PALM / Obs-GWAS)
##   E  split-sample analysis (100 proteins): replication of the discovery-half hits in
##      the validation half (P < 0.05 with the same sign; P < 5e-8), SNPs and primary
##      associations, Syn-PALM and Obs-GWAS
##   F  split-sample analysis: genome-wide significant SNPs per protein in the discovery
##      half, Syn-PALM against Obs-GWAS
## Inputs: results/table_results_cf_noPRS_noPC/external_validation_loci.csv
##         (external_validation.R) and results/split_validation_cf_noPRS_noPC_summary.csv
##         (summarize_split.R).
##   Rscript code/post_cf3/figs/fig_validation_cf.R
## Writes results/figs_cf_noPRS_noPC/fig_validation.{pdf,png} (7.2 x 5 in, print size)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })
setwd("<PROJECT_DIR>")
source("simulation_code/rev1/theme_rev1.R")
OUT <- "results/figs_cf_noPRS_noPC/"
th <- theme_rev1(base_size = 8)
syn_col <- unname(method_cols["Syn-PALM"]); obs_col <- unname(method_cols["Obs-GWAS"])
grp_cols <- c("Obs-GWAS" = obs_col, "Syn-PALM, shared" = "#7FA7D6", "Syn-PALM, additional" = syn_col)
nlp_lab <- function(m) bquote(.(m) ~ -log[10] * italic(P))
lab_box <- function(x, y, hjust, vjust, label)
  annotate("label", x = x, y = y, hjust = hjust, vjust = vjust, size = 2.1, border.colour = NA,
           fill = alpha("white", 0.85), label = label, lineheight = 0.95)

x <- fread("results/table_results_cf_noPRS_noPC/external_validation_loci.csv")
x[, grp := fifelse(method == "ObsG", "Obs-GWAS", fifelse(only, "Syn-PALM, additional", "Syn-PALM, shared"))]
x[, grp := factor(grp, levels = rev(names(grp_cols)))]

## ---- A: pooled shares (horizontal bars) -------------------------------------------
pa <- rbind(x[catalog_records > 0, .(ok = sum(in_catalog), n = .N), by = grp][, what := "Reported in the GWAS Catalog"],
            x[!is.na(sun_sig), .(ok = sum(sun_sig), n = .N), by = grp][, what := "P < 5e-8 in Sun et al."])
pa[, `:=`(share = ok / n, what = factor(what, levels = c("Reported in the GWAS Catalog", "P < 5e-8 in Sun et al.")))]
pa[, `:=`(lab = sprintf("%.1f%%", 100 * share), lab_n = sprintf("%s/%s", comma(ok), comma(n)), inside = share > 0.45)]
cat("A:\n"); print(pa[, .(what, grp, ok, n, share)])
A <- ggplot(pa, aes(share, grp, fill = grp)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = lab), hjust = -0.08, size = 2.1, fontface = "bold") +
  geom_text(data = pa[inside == TRUE], aes(x = 0.02, label = lab_n), hjust = 0, size = 1.9, colour = "white") +
  geom_text(data = pa[inside == FALSE], aes(x = share + 0.25, label = lab_n), hjust = 0, size = 1.9, colour = "grey30") +
  facet_wrap(~ what, ncol = 1) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_x_continuous(labels = percent, limits = c(0, 1.2), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
  scale_y_discrete(labels = c("Obs-GWAS" = "Obs-GWAS", "Syn-PALM, shared" = "Syn-PALM shared",
                              "Syn-PALM, additional" = "Syn-PALM additional")) +
  labs(x = "Primary associations", y = NULL) + th +
  theme(panel.grid.major.y = element_blank())

## ---- B-D: additional Syn-PALM associations ------------------------------------------
so <- x[method == "SynG" & only == TRUE]
so[, Sun := factor(fifelse(!is.na(sun_sig) & sun_sig, "yes", "no"), levels = c("yes", "no"))]
cat(sprintf("additional: %d; Obs P<1e-3 %.1f%%, <1e-5 %.1f%%; median Obs P %.2g; same sign %.1f%%; median beta ratio %.3f; median SE ratio %.3f; Sun sig %.1f%% (of %d with data)\n",
            nrow(so), 100 * mean(so$ObsG_nlp >= 3), 100 * mean(so$ObsG_nlp >= 5), 10^-median(so$ObsG_nlp),
            100 * mean(sign(so$SynSurrG_beta) == sign(so$ObsG_beta)), median(so$SynSurrG_beta / so$ObsG_beta),
            median(so$SynSurrG_se / so$ObsG_se), 100 * mean(so$sun_sig, na.rm = TRUE), sum(!is.na(so$sun_sig))))
thr <- -log10(5e-8)
B <- ggplot(so[order(Sun == "yes")], aes(SynSurrG_nlp, ObsG_nlp, colour = Sun)) +
  geom_hline(yintercept = c(3, 5), linetype = "dotted", colour = "grey35", linewidth = 0.35) +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "red3", linewidth = 0.4) +
  geom_point(size = 0.45, alpha = 0.55, stroke = 0) +
  annotate("text", x = 30, y = c(3, 5), label = c("P == 10^-3", "P == 10^-5"), parse = TRUE,
           hjust = 1, vjust = -0.3, size = 2, colour = "grey30") +
  annotate("text", x = 30, y = thr, label = "P == 5 %*% 10^-8", parse = TRUE, hjust = 1, vjust = -0.3, size = 2, colour = "red3") +
  lab_box(30, 0.15, 1, 0, sprintf("n = %s\nObs-GWAS P < 1e-3: %.1f%%\nObs-GWAS P < 1e-5: %.1f%%",
                                  comma(nrow(so)), 100 * mean(so$ObsG_nlp >= 3), 100 * mean(so$ObsG_nlp >= 5))) +
  scale_x_continuous(trans = "log10", breaks = c(8, 10, 15, 20, 30), limits = c(thr, 30)) +
  scale_y_continuous(limits = c(0, 10.6), breaks = seq(0, 8, 2), expand = c(0, 0)) +
  scale_colour_manual(values = c(yes = syn_col, no = "#B0B7C3"), name = NULL,
                      labels = c(yes = expression(italic(P) < 5 %*% 10^-8 ~ "in Sun et al."), no = "otherwise"),
                      guide = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(x = nlp_lab("Syn-PALM"), y = nlp_lab("Obs-GWAS")) + th +
  theme(legend.position = "inside", legend.position.inside = c(0.01, 0.995), legend.justification = c(0, 1),
        legend.direction = "vertical", legend.background = element_blank(), legend.spacing.y = unit(0, "pt"),
        legend.key.size = unit(0.5, "lines"), legend.text = element_text(size = rel(0.8)))

lim <- quantile(abs(c(so$SynSurrG_beta, so$ObsG_beta)), 0.998)
C <- ggplot(so, aes(SynSurrG_beta, ObsG_beta)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "red3", linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 0.45, alpha = 0.35, stroke = 0, colour = syn_col) +
  lab_box(lim, -lim, 1, 0, sprintf("same sign: %.0f%%\nmedian ratio: %.2f",
                                   100 * mean(sign(so$SynSurrG_beta) == sign(so$ObsG_beta)), median(so$SynSurrG_beta / so$ObsG_beta))) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression("Syn-PALM " * hat(beta)[G]), y = expression("Obs-GWAS " * hat(beta)[G])) + th

ser <- so$SynSurrG_se / so$ObsG_se
D <- ggplot(data.frame(r = ser), aes(r)) +
  geom_histogram(binwidth = 0.0125, boundary = 1, fill = syn_col, colour = "white", linewidth = 0.1) +
  geom_vline(xintercept = 1, colour = "grey30", linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = median(ser), colour = "red3", linewidth = 0.5) +
  annotate("text", x = median(ser), y = Inf, vjust = 1.4, hjust = 1.08, size = 2.1, colour = "red3",
           label = sprintf("median %.3f", median(ser))) +
  scale_x_continuous(limits = quantile(ser, c(0.002, 0.998)) + c(-0.01, 0.01), breaks = seq(0.6, 1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "SE ratio, Syn-PALM / Obs-GWAS", y = "Additional associations") + th

## ---- E-F: split-sample analysis ---------------------------------------------------------
sp <- fread("results/split_validation_cf_noPRS_noPC_summary.csv")
sp[, Method := factor(fifelse(method == "SynSurrG", "Syn-PALM", "Obs-GWAS"), levels = c("Syn-PALM", "Obs-GWAS"))]
pe <- rbind(sp[, .(level = "GWS SNPs", crit = "nom", k = sum(snp_nom), n = sum(snp_disc)), by = Method],
            sp[, .(level = "GWS SNPs", crit = "gw", k = sum(snp_overlap), n = sum(snp_disc)), by = Method],
            sp[, .(level = "Primary associations", crit = "nom", k = sum(locus_nom), n = sum(locus_disc)), by = Method],
            sp[, .(level = "Primary associations", crit = "gw", k = sum(locus_overlap), n = sum(locus_disc)), by = Method])
pe[, `:=`(share = k / n, crit = factor(crit, levels = c("nom", "gw")))]
pe[, level := sprintf("%s\n%s vs %s", sub("Primary associations", "Primary assoc.", level), comma(n[Method == "Syn-PALM"][1]), comma(n[Method == "Obs-GWAS"][1])), by = level]
cat("E:\n"); print(pe)
E <- ggplot(pe, aes(crit, share, fill = Method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * share)), position = position_dodge(0.8), angle = 90,
            hjust = 1.1, size = 2, colour = "white") +
  facet_wrap(~ level) +
  scale_fill_manual(values = method_cols[c("Syn-PALM", "Obs-GWAS")], name = NULL) +
  scale_x_discrete(labels = c(nom = "P < 0.05,\nsame sign", gw = "P < 5e-8")) +
  scale_y_continuous(labels = percent, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Replication in the validation half", y = "Discovery hits replicated") + th +
  theme(legend.position = "bottom", legend.key.size = unit(0.6, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "pt"))

pf <- dcast(sp, protein ~ Method, value.var = "snp_disc")
cat(sprintf("F: SNPs discovered Syn %d, Obs %d; proteins with more Syn %d of %d\n",
            sum(pf$`Syn-PALM`), sum(pf$`Obs-GWAS`), sum(pf$`Syn-PALM` > pf$`Obs-GWAS`), nrow(pf)))
fmax <- max(c(pf$`Syn-PALM`, pf$`Obs-GWAS`)) * 1.3
F <- ggplot(pf, aes(`Obs-GWAS` + 1, `Syn-PALM` + 1)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = log10(2), colour = "grey60", linetype = "dotted", linewidth = 0.4) +
  geom_point(size = 1.1, alpha = 0.75, colour = syn_col) +
  lab_box(fmax, 1, 1, 0, sprintf("%d of %d proteins\nabove y = x", sum(pf$`Syn-PALM` > pf$`Obs-GWAS`), nrow(pf))) +
  scale_x_log10(limits = c(1, fmax)) + scale_y_log10(limits = c(1, fmax)) + coord_equal() +
  labs(x = "Obs-GWAS GWS SNPs + 1", y = "Syn-PALM GWS SNPs + 1") + th

## ---- layout -------------------------------------------------------------------------
tag <- theme(plot.tag = element_text(face = "bold", size = 10))
## each row is wrapped so that its axes are not aligned with those of the other row
## (A has long y labels); tags are set per panel
tg <- function(p, l) p + labs(tag = l) + tag
fig <- wrap_elements(full = tg(A, "A") + tg(B, "B") + tg(C, "C") + plot_layout(widths = c(1.2, 1.25, 0.9))) /
       wrap_elements(full = tg(D, "D") + tg(E, "E") + tg(F, "F") + plot_layout(widths = c(0.95, 1.45, 0.85)))
ggsave(paste0(OUT, "fig_validation.pdf"), fig, width = 7.2, height = 5.0, device = cairo_pdf)
ggsave(paste0(OUT, "fig_validation.png"), fig, width = 7.2, height = 5.0, dpi = 300)
cat("wrote", paste0(OUT, "fig_validation.pdf"), "\n")

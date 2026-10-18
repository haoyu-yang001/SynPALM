## Fig 4F and Supplementary Figure obssea (revision): SNPsea tissue enrichment of the
## lead variants of Syn-PALM and Obs-GWAS (code/sh_run/run_noPRS_noPC/run_snpsea.sh),
## as in the submitted panels: the 30 Gene Atlas tissues with the smallest P values,
## ordered by UPGMA clustering of the Pearson correlations of their expression
## profiles (heatmap), with -log10 P bars and the Bonferroni line (0.05 / 79).
##   Rscript code/post_cf3/figs/fig4f_snpsea.R
## Writes results/figs_cf_noPRS_noPC/snpsea_{SynG,ObsG,SynOnly}.pdf and
## results/table_results_cf_noPRS_noPC/snpsea_condition_pvalues.csv
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
setwd("<PROJECT_DIR>")
S <- "results/table_results_cf_noPRS_noPC/snpsea/"
gct <- fread(cmd = paste("zcat", path.expand("~/tools/snpsea_data/GeneAtlas2004.gct.gz")), skip = 2)
E <- as.matrix(gct[, -(1:2)])
lab <- c(SynG = "Syn-PALM", ObsG = "Obs-GWAS", SynOnly = "Syn-PALM only")
all <- rbindlist(lapply(names(lab), function(s) {
  f <- paste0(S, s, "/condition_pvalues.txt"); if (!file.exists(f)) return(NULL)
  cbind(set = s, fread(f))
}))
fwrite(all, "results/table_results_cf_noPRS_noPC/snpsea_condition_pvalues.csv")
bonf <- -log10(0.05 / ncol(E))
for (s in unique(all$set)) {
  x <- all[set == s][order(pvalue)][1:30]
  C <- cor(E[, x$condition])
  o <- hclust(as.dist(1 - C), method = "average")$order
  lev <- x$condition[o]
  hm <- as.data.table(as.table(C))[, `:=`(V1 = factor(V1, lev), V2 = factor(V2, lev))]
  hm <- hm[as.integer(V1) >= as.integer(V2)]
  p1 <- ggplot(hm, aes(V2, V1, fill = N)) + geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Pearson's r") +
    theme_void() + theme(legend.position = "bottom", legend.key.height = unit(0.25, "cm"))
  x[, condition := factor(condition, lev)]
  x[, sig := -log10(pvalue) > bonf]
  p2 <- ggplot(x, aes(-log10(pvalue), condition, fill = sig)) + geom_col(width = 0.7) +
    geom_vline(xintercept = bonf, colour = "black") +
    scale_fill_manual(values = c(`TRUE` = "#3399FF", `FALSE` = "grey65"), guide = "none") +
    scale_y_discrete(position = "right", labels = function(v) ifelse(x$sig[match(v, x$condition)], paste("*", v), v)) +
    labs(x = expression(-log[10] * italic(P)), y = NULL, title = lab[[s]]) +
    theme_classic() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5, face = "bold"))
  p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(0.8, 1), align = "h", axis = "tb")
  ggsave(sprintf("results/figs_cf_noPRS_noPC/snpsea_%s.pdf", s), p, width = 7, height = 5.5)
  cat(s, ": tissues with P < 0.05/79:", paste(sprintf("%s (%.1e)", x[sig == TRUE, condition], x[sig == TRUE, pvalue]), collapse = "; "), "\n")
}

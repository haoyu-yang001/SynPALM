## Fig 4C (revision): primary associations (P < 5e-8) of Obs-GWAS and Syn-PALM for the
## 20 proteins with the highest out-of-fold prediction correlation (in parentheses),
## as in the submitted panel; counts from aggregate_all.R (summary_per_protein.csv).
##   Rscript code/post_cf3/figs/fig4c_top20_cf.R
## Writes results/figs_cf_noPRS_noPC/fig4c_top20.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setwd("<PROJECT_DIR>")
per <- fread("results/table_results_cf_noPRS_noPC/summary_per_protein.csv")
top <- per[order(-rho_oof)][1:20]
top[, lab := sprintf("%s (%.3f)", tolower(protein), rho_oof)]
d <- melt(top[, .(lab, `Obs-GWAS` = `ObsG_loci_5e-08`, `Syn-PALM` = `SynG_loci_5e-08`)], id.vars = "lab",
          variable.name = "Method", value.name = "n")
d[, lab := factor(lab, levels = top$lab)]
p <- ggplot(d, aes(lab, n, fill = Method)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c(`Obs-GWAS` = "#C0392B", `Syn-PALM` = "#0073C2")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Proteins", y = "GWS associations", fill = "Method") + theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7), legend.position = "top")
ggsave("results/figs_cf_noPRS_noPC/fig4c_top20.pdf", p, width = 7.5, height = 4)

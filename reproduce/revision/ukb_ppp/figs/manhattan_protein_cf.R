## Revision copy of code/gwas_plot_for_protein.R (Fig 5C): Manhattan plots of
## Syn-PALM and Obs-GWAS for one protein from the no-PRS / no-PC rerun, with the
## lead variants of the genome-wide significant loci (P < 5e-8) labelled by
## their nearest gene (TxDb hg38).
##   Rscript code/post_cf3/figs/manhattan_protein_cf.R [protein=gdf15]
## Writes results/figs_cf_noPRS_noPC/Manhattan_{obs,syn}_gwas_<protein>.{pdf,png}
suppressPackageStartupMessages({ library(dplyr); library(scales) })
setwd("<PROJECT_DIR>")
source("code/gwas_plot_function.R")
protein_name <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "gdf15"
OUT <- "results/figs_cf_noPRS_noPC/"; dir.create(OUT, showWarnings = FALSE)
thr <- -log10(5e-8)
pv <- readRDS(sprintf("results/web_SynGWAS_summary_results_cf_noPRS_noPC/%s_negative_log10_pvals_v2geno_crossfit.rds", protein_name))
L  <- readRDS(sprintf("results/table_results_cf_noPRS_noPC/%s_loci.rds", protein_name))
for (m in c("ObsG", "SynG")) {
  col <- if (m == "ObsG") "ObsG_negative_log10_pvalues" else "SynSurrG_negative_log10_pvalues"
  d <- data.frame(rsid = pv$rsid, chr = pv$chr, logP = pv[[col]], pos = pv$POS38)
  lead <- L$by_threshold[["5e-08"]][[m]]$loci
  hl <- if (nrow(lead)) {
    ng <- get_nearest_gene_from_pos(data.frame(rsid = lead$rsid, chr = lead$chr, pos = lead$POS38), build = "hg38")
    data.frame(rsid = ng$rsid, gene = ng$nearest_gene)
  } else data.frame(rsid = character(0), gene = character(0))
  print(hl)
  p <- plot_manhattan_refined(data = d, p_col = "logP", title_name = if (m == "ObsG") "Obs-GWAS" else "Syn-PALM",
                              highlight_df = hl, thre = thr)
  ## fixed breaks: the automatic log breaks put 0 and 0.x on top of each other (Reviewer 3, minor 5)
  p <- suppressMessages(p + scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                                               breaks = c(0, 1, 10, 100, 1000), labels = label_comma(accuracy = 1),
                                               expand = expansion(mult = c(0, 0.15))))
  ## revision: gene labels lifted above the significance line (they sat on the dashed line),
  ## chromosomes 18, 20 and 21 not labelled (the labels of the short chromosomes ran
  ## together; 19 is kept because GDF15 is on it), larger text for the composite figure
  is_lab <- vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1))
  lab_df <- p$layers[[which(is_lab)]]$data
  p$layers <- p$layers[!is_lab]
  ctr <- p$data %>% group_by(chr) %>% summarise(center = (max(BPcum) + min(BPcum)) / 2)
  p <- suppressMessages(p +
    geom_text_repel(data = lab_df, aes(label = gene), size = 3.6, fontface = "italic", color = "black",
                    nudge_y = 0.35, direction = "both", min.segment.length = 0, box.padding = 0.4,
                    point.padding = 0.15, force = 8, segment.color = "grey40", segment.size = 0.35,
                    max.overlaps = Inf, seed = 42) +
    scale_x_continuous(breaks = ctr$center, labels = ifelse(ctr$chr %in% c(18, 20, 21), "", ctr$chr),
                       expand = expansion(mult = 0.01)) +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12), plot.title = element_text(size = 14)))
  tag <- if (m == "ObsG") "obs" else "syn"
  for (ext in c("pdf", "png"))
    ggsave(sprintf("%sManhattan_%s_gwas_%s.%s", OUT, tag, protein_name, ext), p, width = 6.5, height = 3.4, dpi = 300)
}

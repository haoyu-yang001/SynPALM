## Fig 5B (revision): Circos plot of the SuSiE fine-mapping results of Syn-PALM
## (code/post_cf3/finemap_susie.R) for the 25 proteins with the highest
## out-of-fold prediction correlation. Revision of code/finemap_results/finemap_plot_code.R:
## the credible sets come from susie_rss instead of single-variant ABFs, the
## inner track shows the Syn-PALM -log10 P of the variants with PIP >= 0.01
## (the old code plotted the Obs-GWAS P values), the links go from the gene of
## the protein to the top variant of every 95% credible set with PIP >= 0.95, and
## only trans signals (more than 1 Mb from the gene or on another chromosome) are linked;
## cis signals sit next to the gene label.
## As in the submitted figure, one variant is labelled and linked per protein and
## chromosome: the credible-set top variant with the highest PIP (if >= 0.95).
## The inner track is truncated at -log10 P = 50.
##   Rscript code/post_cf3/figs/fig5b_circos_susie.R [top_n=25]
## Writes results/figs_cf_noPRS_noPC/fig5b_circos_susie.{pdf,png} and
## results/table_results_cf_noPRS_noPC/fig5b_circos_susie_links.csv
suppressPackageStartupMessages({ library(circlize); library(data.table) })
setwd("<PROJECT_DIR>")
TOP_N <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 25L
SNP_LABEL_PIP <- 0.95; PIP_LINK <- 0.95
SNP_LABEL_CEX <- 0.38; PROT_LABEL_CEX <- 0.45
TRANS_W <- 2e5; CIS_W <- 4e5
OUT <- "results/figs_cf_noPRS_noPC/"

rho <- as.data.table(readRDS("results/rho_oof_cf_noPRS_noPC.rds"))[order(-rho_oof)]
sel <- rho$protein[seq_len(TOP_N)]
gmap <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
gmap[, key := gsub("-", "_", id)]

S <- lapply(sel, function(p) readRDS(sprintf("results/table_results_cf_noPRS_noPC/susie/%s.rds", p))$SynG)
names(S) <- sel
var <- rbindlist(lapply(S, `[[`, "var"), fill = TRUE)
cs  <- rbindlist(lapply(S, `[[`, "cs"), fill = TRUE)
var[, `:=`(chr_c = paste0("chr", chr), nlp = -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10), in_cs = !is.na(cs))]

pal <- colorRampPalette(c("#FF6B9D", "#FF8E53", "#FFD93D", "#6BCB77", "#4ECDC4", "#45B7D1", "#A29BFE", "#FD79A8",
                          "#00CEC9", "#FDCB6E", "#E17055", "#74B9FF", "#55EFC4", "#B8E994", "#E84393", "#00B894",
                          "#0984E3", "#6C5CE7", "#F8C471", "#82E0AA"))(TOP_N)
prot_pal <- setNames(pal, sel)

## gene positions of the selected proteins (autosomes)
gp <- gmap[key %in% sel][!duplicated(key)][as.character(chr) %in% as.character(1:22)]
gp[, `:=`(chr_c = paste0("chr", chr), color = prot_pal[key], label = toupper(gsub("_", "/", key)))]

## top variant of each credible set
top <- cs[, .(protein, chr_c = paste0("chr", chr), pos = top_pos, rsid = top_rsid, pip = top_pip)]
top[, color := prot_pal[protein]]
best <- top[order(-pip)][!duplicated(paste(protein, chr_c))]
lab <- best[pip >= SNP_LABEL_PIP][order(chr_c, pos)][!duplicated(rsid)]
links <- merge(best[pip >= PIP_LINK], gp[, .(protein = key, gene_chr = chr_c, gene_pos = pos)], by = "protein")
links[, is_cis := gene_chr == chr_c & abs(gene_pos - pos) < 1e6]
links <- links[is_cis == FALSE]      # cis signals (within 1 Mb of the gene) cannot be drawn as arcs; they sit next to the gene label
fwrite(links, "results/table_results_cf_noPRS_noPC/fig5b_circos_susie_links.csv")
cat(sprintf("proteins %d (rho %.3f-%.3f); credible sets %d; top variants with PIP >= %.2f: %d; links %d (cis %d, trans %d)\n",
            TOP_N, min(rho$rho_oof[seq_len(TOP_N)]), max(rho$rho_oof[seq_len(TOP_N)]), nrow(cs), PIP_LINK,
            top[pip >= PIP_LINK, .N], nrow(links), links[is_cis == TRUE, .N], links[is_cis == FALSE, .N]))
print(cs[, .(credible_sets = .N, chromosomes = uniqueN(chr), pip95 = sum(top_pip >= 0.95)), by = protein][order(-credible_sets)])

darken <- function(col, f) { v <- col2rgb(col) / 255; rgb(v[1, ] * f, v[2, ] * f, v[3, ] * f) }
draw <- function() {
  circos.par(start.degree = 90, gap.degree = c(rep(0.8, 21), 3), track.margin = c(0.003, 0.003),
             cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE, circle.margin = 0.45)
  circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", 1:22),
                                plotType = c("ideogram", "labels"), ideogram.height = 0.03, labels.cex = 0.55)
  if (nrow(lab))
    circos.genomicLabels(data.frame(chr = lab$chr_c, start = lab$pos, end = lab$pos + 1, label = lab$rsid),
                         labels.column = 4, side = "outside", col = lab$color, line_col = lab$color,
                         cex = SNP_LABEL_CEX, font = 1, line_lwd = 0.6, padding = 0.5)
  ymax <- 52
  circos.track(track.height = 0.09, bg.col = "#F8F8F8", bg.border = "#DDDDDD", ylim = c(0, ymax),
    panel.fun = function(region, value, ...) {
      sub <- var[chr_c == get.cell.meta.data("sector.index")]
      if (nrow(sub)) for (i in seq_len(nrow(sub)))
        circos.lines(rep(sub$POS38[i], 2), c(0, min(sub$nlp[i], ymax)),
                     col = if (sub$in_cs[i]) prot_pal[sub$protein[i]] else if (sub$z[i] > 0) "#FFCCCC" else "#CCCCFF",
                     lwd = if (sub$in_cs[i]) 1.5 else 0.4)
      circos.lines(get.cell.meta.data("xlim"), rep(-log10(5e-8), 2), col = "red", lty = 2, lwd = 0.5)
    })
  if (nrow(gp))
    circos.genomicLabels(data.frame(chr = gp$chr_c, start = gp$pos, end = gp$pos + 1, label = gp$label),
                         labels.column = 4, side = "inside", col = gp$color, line_col = gp$color,
                         cex = PROT_LABEL_CEX, font = 2, line_lwd = 0.8, padding = 0.5)
  for (i in seq_len(nrow(links))) {
    l <- links[i]
    tryCatch(circos.link(l$gene_chr, l$gene_pos, l$chr_c, l$pos,
                         col = adjustcolor(darken(l$color, if (l$is_cis) 0.5 else 0.85), alpha.f = if (l$is_cis) 1 else 0.75),
                         lwd = if (l$is_cis) 2.5 else 1, h.ratio = if (l$is_cis) 0.35 else 0.6), error = function(e) NULL)
  }
  legend(x = 0.72, y = -0.68, bty = "n", cex = 0.48,
         legend = c("In a 95% credible set", "Other variant with PIP >= 0.01", expression(P == 5 %*% 10^-8),
                    "trans-pQTL (gene -> variant)"),
         col = c("black", "#CCCCCC", "red", "#888888"), pch = c(18, 16, NA, NA),
         lty = c(NA, NA, 2, 1), lwd = c(NA, NA, 1, 1), title = expression("Inner track: Syn-PALM " * -log[10] * P * " (truncated at 50)"))
  circos.clear()
}
cairo_pdf(paste0(OUT, "fig5b_circos_susie.pdf"), width = 12, height = 12); draw(); invisible(dev.off())
png(paste0(OUT, "fig5b_circos_susie.png"), width = 12, height = 12, units = "in", res = 200, type = "cairo"); draw(); invisible(dev.off())
cat("wrote", paste0(OUT, "fig5b_circos_susie.pdf"), "\n")

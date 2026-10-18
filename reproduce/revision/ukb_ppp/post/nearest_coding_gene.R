## Nearest protein-coding gene (UCSC knownGene hg38, Entrez gene type
## "protein-coding" in org.Hs.eg.db) of each variant; distance 0 inside a gene.
##   source("code/post_cf3/nearest_coding_gene.R"); nearest_coding_gene(rsid, chr, pos)
nearest_coding_gene <- function(rsid, chr, pos) {
  suppressPackageStartupMessages({ library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db) })
  g <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  ann <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = names(g), keytype = "ENTREZID", columns = c("SYMBOL", "GENETYPE")))
  ann <- ann[!duplicated(ann$ENTREZID), ]
  keep <- ann$ENTREZID[ann$GENETYPE %in% "protein-coding"]
  g <- g[names(g) %in% keep]
  g <- GenomeInfoDb::keepSeqlevels(g, paste0("chr", 1:22), pruning.mode = "coarse")
  s <- GRanges(paste0("chr", sub("^chr", "", as.character(chr))), IRanges(as.numeric(pos), width = 1))
  i <- nearest(s, g, ignore.strand = TRUE)
  out <- data.frame(rsid = rsid, nearest_gene = NA_character_, distance = NA_integer_)
  ok <- !is.na(i)
  out$nearest_gene[ok] <- ann$SYMBOL[match(names(g)[i[ok]], ann$ENTREZID)]
  out$distance[ok] <- distance(s[ok], g[i[ok]], ignore.strand = TRUE)
  out
}

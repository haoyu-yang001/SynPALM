## Split-sample replication summary for split_validation_crossfit.R output.
##   Rscript code/post_cf3/summarize_split.R
## Per protein and method (SynSurrG = Syn-PALM, ObsG = observed-only LMM):
##   snp_overlap : discovery hits (p < 5e-8) also p < 5e-8 in validation
##                 (the submitted manuscript's definition)
##   snp_nom     : discovery hits with p < 0.05 and the same sign in validation
##   locus_*     : the same after clumping discovery hits into +-1 Mb loci
##                 (lead = smallest p); a locus replicates if its lead SNP does
## Autosomes only; MHC (chr6:25.5-34 Mb, GRCh38) is left in; SNPs without POS38 use POS37.
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
setwd("<PROJECT_DIR>")
D   <- "results/split_validation_cf_noPRS_noPC/"
ann <- readRDS("data/snp_annot_v2geno.rds")
THR <- -log10(5e-8)

clump_leads <- function(chr, pos, lp) {
  o <- order(-lp); keep <- integer(0); used <- logical(length(o))
  for (i in o) {
    if (used[i]) next
    keep <- c(keep, i)
    used[chr == chr[i] & abs(pos - pos[i]) <= 1e6] <- TRUE
  }
  keep
}

files <- list.files(D, pattern = "_split_scan\\.rds$", full.names = TRUE)
rows <- list()
for (f in files) {
  x <- readRDS(f); p <- sub("_split_scan\\.rds$", "", basename(f))
  a <- ann[match(rownames(x), ann$rsid), ]
  pos <- ifelse(is.na(a$POS38), a$POS37, a$POS38)
  auto <- a$chr %in% as.character(1:22)
  for (m in c("SynSurrG", "ObsG")) {
    ld <- x[, paste0(m, "_discovery_negative_log10_pvalues")]
    lv <- x[, paste0(m, "_validation_negative_log10_pvalues")]
    bd <- x[, paste0(m, "_discovery_hat_beta")]; bv <- x[, paste0(m, "_validation_hat_beta")]
    ld[!auto] <- NA; lv[!auto] <- NA
    hit <- which(ld >= THR)
    rep_gw  <- lv[hit] >= THR
    rep_nom <- lv[hit] >= -log10(0.05) & sign(bd[hit]) == sign(bv[hit])
    lead <- if (length(hit)) hit[clump_leads(a$chr[hit], pos[hit], ld[hit])] else integer(0)
    rows[[length(rows) + 1]] <- data.frame(
      protein = p, method = m,
      rho_disc = attr(x, "rho")[["discovery"]], rho_val = attr(x, "rho")[["validation"]],
      n_obs_disc = attr(x, "n_obs")[["discovery"]],
      snp_disc = length(hit), snp_val = sum(lv >= THR, na.rm = TRUE),
      snp_overlap = sum(rep_gw, na.rm = TRUE), snp_nom = sum(rep_nom, na.rm = TRUE),
      locus_disc = length(lead),
      locus_overlap = sum(lv[lead] >= THR, na.rm = TRUE),
      locus_nom = sum(lv[lead] >= -log10(0.05) & sign(bd[lead]) == sign(bv[lead]), na.rm = TRUE))
  }
}
tab <- bind_rows(rows)
saveRDS(tab, "results/split_validation_cf_noPRS_noPC_summary.rds")
fwrite(tab, "results/split_validation_cf_noPRS_noPC_summary.csv")

cat("proteins:", length(files), "\n")
for (m in c("SynSurrG", "ObsG")) {
  t <- filter(tab, method == m, snp_disc > 0)
  cat(sprintf("\n%s (%d proteins with >=1 discovery hit)\n", m, nrow(t)))
  cat(sprintf("  total discovery SNP hits %d, loci %d\n", sum(t$snp_disc), sum(t$locus_disc)))
  q <- function(v) paste(sprintf("%.1f%%", 100 * quantile(v, c(.25, .5, .75))), collapse = " / ")
  cat("  per-protein SNP replication p<5e-8   Q1/median/Q3:", q(t$snp_overlap / t$snp_disc),
      sprintf(" mean %.1f%%; pooled %.1f%%\n", 100 * mean(t$snp_overlap / t$snp_disc), 100 * sum(t$snp_overlap) / sum(t$snp_disc)))
  cat("  per-protein SNP replication p<0.05   Q1/median/Q3:", q(t$snp_nom / t$snp_disc),
      sprintf(" pooled %.1f%%\n", 100 * sum(t$snp_nom) / sum(t$snp_disc)))
  cat(sprintf("  locus replication pooled: p<5e-8 %.1f%%, p<0.05 same sign %.1f%%\n",
              100 * sum(t$locus_overlap) / sum(t$locus_disc), 100 * sum(t$locus_nom) / sum(t$locus_disc)))
}

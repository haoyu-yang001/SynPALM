## Reply to Reviewer 3, Comment 3 (revision 2026-10-20): replication of the additional
## Syn-PALM associations (no Obs-GWAS primary association within +-1 Mb) against
## comparators of similar strength.
## (A) External (external_validation_loci.csv, all 2,916 proteins): share reported in the
##     GWAS Catalog and share with P < 5e-8 in Sun et al., for Syn-PALM loci (shared /
##     additional) in bins of their Syn-PALM -log10 P, and for Obs-GWAS loci in bins of
##     their Obs-GWAS -log10 P; additional loci also by their Obs-GWAS P.
## (B) Internal (split-sample scans of the 100 proteins, results/split_validation_cf_noPRS_noPC/):
##     discovery-half Syn-PALM loci (greedy +-1 Mb clumping, as summarize_split.R) split into
##     additional (no Obs-GWAS discovery hit within +-1 Mb) and shared, with replication in
##     the validation half (P < 0.05 and same sign; P < 5e-8), by bins of the discovery P;
##     Obs-GWAS discovery loci by the same bins as the comparator.
##   Rscript code/post_cf3/additional_loci_replication.R
## Writes results/table_results_cf_noPRS_noPC/additional_loci_replication_{external,split}.csv
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
OUT <- "results/table_results_cf_noPRS_noPC/"
BR <- c(-log10(5e-8), 8, 9, 10, 15, Inf); BL <- c("5e-8 to 1e-8", "1e-8 to 1e-9", "1e-9 to 1e-10", "1e-10 to 1e-15", "< 1e-15")
x <- fread(paste0(OUT, "external_validation_loci.csv"))
x[, group := fifelse(method == "ObsG", "Obs-GWAS", fifelse(only, "Syn-PALM additional", "Syn-PALM shared"))]
x[, own_nlp := fifelse(method == "ObsG", ObsG_nlp, SynSurrG_nlp)]
x[, bin := cut(own_nlp, BR, BL, right = FALSE)]
ext <- x[, .(loci = .N, catalog = mean(in_catalog[catalog_records > 0]), n_cat = sum(catalog_records > 0),
             sun = mean(sun_sig, na.rm = TRUE)), by = .(group, bin)][order(group, bin)]
cat("(A) by own P bin\n"); print(ext, digits = 3)
add <- x[group == "Syn-PALM additional"]
add[, obs_bin := cut(ObsG_nlp, c(-Inf, 5, 6, 7, -log10(5e-8)), c("> 1e-5", "1e-5 to 1e-6", "1e-6 to 1e-7", "1e-7 to 5e-8"), right = FALSE)]
ext2 <- add[, .(loci = .N, catalog = mean(in_catalog[catalog_records > 0]), sun = mean(sun_sig, na.rm = TRUE)), by = obs_bin][order(obs_bin)]
cat("(A) additional loci by Obs-GWAS P\n"); print(ext2, digits = 3)
fwrite(rbind(ext[, .(table = "own_bin", group, bin = as.character(bin), loci, catalog, sun)],
             ext2[, .(table = "additional_by_obsP", group = "Syn-PALM additional", bin = as.character(obs_bin), loci, catalog, sun)]),
       paste0(OUT, "additional_loci_replication_external.csv"))
## (B) split sample
ann <- readRDS("data/snp_annot_v2geno.rds")
THR <- -log10(5e-8)
clump_leads <- function(chr, pos, lp) {
  o <- order(-lp); keep <- integer(0); used <- logical(length(o))
  for (i in o) { if (used[i]) next; keep <- c(keep, i); used[chr == chr[i] & abs(pos - pos[i]) <= 1e6] <- TRUE }
  keep
}
rows <- list()
for (f in list.files("results/split_validation_cf_noPRS_noPC", "_split_scan\\.rds$", full.names = TRUE)) {
  s <- readRDS(f); a <- ann[match(rownames(s), ann$rsid), ]
  pos <- ifelse(is.na(a$POS38), a$POS37, a$POS38); auto <- a$chr %in% as.character(1:22)
  get <- function(m, h, k) { v <- s[, sprintf("%s_%s_%s", m, h, k)]; v[!auto] <- NA; v }
  hitO <- which(get("ObsG", "discovery", "negative_log10_pvalues") >= THR)
  for (m in c("SynSurrG", "ObsG")) {
    ld <- get(m, "discovery", "negative_log10_pvalues"); hit <- which(ld >= THR)
    if (!length(hit)) next
    lead <- hit[clump_leads(a$chr[hit], pos[hit], ld[hit])]
    lv <- get(m, "validation", "negative_log10_pvalues")
    bd <- get(m, "discovery", "hat_beta"); bv <- get(m, "validation", "hat_beta")
    addl <- vapply(lead, function(i) !any(a$chr[hitO] == a$chr[i] & abs(pos[hitO] - pos[i]) <= 1e6), logical(1))
    rows[[length(rows) + 1]] <- data.table(protein = sub("_split_scan\\.rds$", "", basename(f)), method = m,
      additional = if (m == "SynSurrG") addl else NA, disc_nlp = ld[lead],
      rep_nom = lv[lead] >= -log10(0.05) & sign(bd[lead]) == sign(bv[lead]), rep_gw = lv[lead] >= THR)
  }
}
sp <- rbindlist(rows)
sp[, group := fifelse(method == "ObsG", "Obs-GWAS", fifelse(additional, "Syn-PALM additional", "Syn-PALM shared"))]
sp[, bin := cut(disc_nlp, BR, BL, right = FALSE)]
tabB <- rbind(sp[, .(bin = "all", loci = .N, rep_nom = mean(rep_nom), rep_gw = mean(rep_gw)), by = group],
              sp[, .(loci = .N, rep_nom = mean(rep_nom), rep_gw = mean(rep_gw)), by = .(group, bin = as.character(bin))])[order(group, bin)]
cat("(B) split sample\n"); print(tabB, digits = 3)
fwrite(tabB, paste0(OUT, "additional_loci_replication_split.csv"))

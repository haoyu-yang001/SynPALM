## Cross-protein summary of the no-PRS / no-PC rerun (reads the <p>_loci.rds
## files written by post_protein.R).
##   Rscript code/post_cf3/aggregate_all.R
## Writes to results/table_results_cf_noPRS_noPC/:
##   summary_per_protein.csv  per protein: rho, n_obs, significant SNPs (all /
##                            non-MHC) and loci for Syn-PALM and Obs-GWAS, 4 thresholds
##   loci_SynG_5e-08.csv ...  lead SNPs with cis/trans (cis: same chromosome as the
##                            protein's gene and within 1 Mb of the gene position)
## and the web-format table (same columns as the submitted
## results/table_for_web/significant_loci_v2geno_all_crossfit_noPRS.csv) to
## results/table_for_web_cf_noPRS_noPC/significant_loci_v2geno_all_cf_noPRS_noPC.csv.
## and its portal-format copy to results/web_SynGWAS_summary_results_cf_noPRS_noPC/.
## Nothing is written to the web portal's directory.
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
D   <- "results/table_results_cf_noPRS_noPC/"
WEB <- "results/table_for_web_cf_noPRS_noPC/"
dir.create(WEB, showWarnings = FALSE)
THR  <- c("5e-08", "1e-08", "1e-09", "1e-10")

files <- list.files(D, pattern = "_loci\\.rds$", full.names = TRUE)
L <- lapply(files, readRDS)
names(L) <- vapply(L, `[[`, "", "protein")
cat("proteins with results:", length(L), "\n")

gmap <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
gmap[, key := gsub("-", "_", id)]
pn <- fread("results/table_results/protein_names.csv", header = TRUE)
setnames(pn, 1:2, c("id", "protein_full_name"))
pn[, key := gsub("-", "_", tolower(sub(":.*", "", id)))]

## ---- per-protein table ------------------------------------------------
per <- rbindlist(lapply(L, function(x) {
  r <- data.table(protein = x$protein, rho_oof = x$rho_oof, n_obs = x$n_obs)
  for (t in THR) {
    b <- x$by_threshold[[t]]
    r[, paste0("SynG_snps_noMHC_", t) := b$SynG$n_snps]
    r[, paste0("ObsG_snps_noMHC_", t) := b$ObsG$n_snps]
    r[, paste0("SynG_loci_", t) := nrow(b$SynG$loci)]
    r[, paste0("ObsG_loci_", t) := nrow(b$ObsG$loci)]
  }
  if (!is.null(x$n_snps_all)) for (i in seq_along(THR)) {
    r[, paste0("SynG_snps_all_", THR[i]) := x$n_snps_all[i, "SynG"]]
    r[, paste0("ObsG_snps_all_", THR[i]) := x$n_snps_all[i, "ObsG"]]
  }
  r
}), fill = TRUE)
setorder(per, -rho_oof)
fwrite(per, paste0(D, "summary_per_protein.csv"))

## ---- lead SNPs with cis / trans ----------------------------------------
cis_trans <- function(loci) {
  g <- gmap[match(gsub("-", "_", loci$protein), key)]
  loci[, `:=`(gene_chr = as.character(g$chr), gene_pos = g$pos)]
  loci[, cis_trans := fifelse(is.na(gene_chr), NA_character_,
                       fifelse(as.character(chr) == gene_chr & abs(POS38 - gene_pos) < 1e6, "cis", "trans"))]
  loci
}
lead <- list()
for (t in THR) for (m in c("SynG", "ObsG")) {
  z <- rbindlist(lapply(L, function(x) {
    lo <- x$by_threshold[[t]][[m]]$loci
    if (nrow(lo)) cbind(protein = x$protein, rho_oof = x$rho_oof, lo) else NULL
  }), fill = TRUE)
  if (nrow(z)) z <- cis_trans(z)
  lead[[paste(m, t)]] <- z
  fwrite(z, paste0(D, "loci_", m, "_", t, ".csv"))
}

## ---- web-format table: union of Syn-PALM and Obs-GWAS leads per protein ----
fmt_p <- function(lp) {
  p <- 10^(-lp); out <- formatC(p, format = "e", digits = 3)
  out[lp > 200] <- "<1e-200"; out
}
web <- rbindlist(lapply(THR, function(t) {
  u <- rbind(lead[[paste("SynG", t)]], lead[[paste("ObsG", t)]], fill = TRUE)
  if (!nrow(u)) return(NULL)
  u <- unique(u, by = c("protein", "rsid"))
  data.table(chr = u$chr, rsid = u$rsid, POS38 = u$POS38, REF = u$REF, ALT = u$ALT,
             EA = u$A1,
             SynGWAS_hat_beta = u$SynSurrG_hat_beta, SynGWAS_sd_hat_beta = sqrt(u$SynSurrG_var_hat_beta),
             SynGWAS_pvalues = fmt_p(u$SynSurrG_negative_log10_pvalues),
             ObsGWAS_hat_beta = u$ObsG_hat_beta, ObsGWAS_sd_hat_beta = sqrt(u$ObsG_var_hat_beta),
             ObsGWAS_pvalues = fmt_p(u$ObsG_negative_log10_pvalues),
             protein_short_name = u$protein, threshold = as.numeric(t),
             predicted_correlation = u$rho_oof,
             protein_full_name = pn$protein_full_name[match(gsub("-", "_", u$protein), pn$key)],
             cis_trans = u$cis_trans)
}))
setorder(web, protein_short_name, -threshold, chr, POS38)
fwrite(web, paste0(WEB, "significant_loci_v2geno_all_cf_noPRS_noPC.csv"))
## the same table in the portal's format (columns, file name and threshold values
## of hyyang_hufeng/SynGWAS_summary_results/significant_loci_v2geno_all_crossfit_noPRS.csv),
## next to the per-protein summary statistics written by post_protein.R
portal <- copy(web)[, c("EA", "cis_trans") := NULL]
portal[, chr := as.integer(chr)]
portal[, threshold := 10^(-as.numeric(as.character(-log10(threshold))))]   # old: exponent parsed from file names
fwrite(portal, "results/web_SynGWAS_summary_results_cf_noPRS_noPC/significant_loci_v2geno_all_crossfit_noPRS.csv")

## ---- headline numbers -------------------------------------------------
s <- function(col) sum(per[[col]], na.rm = TRUE)
cat(sprintf("\nmedian out-of-fold rho: %.3f (n = %d proteins)\n", median(per$rho_oof), nrow(per)))
for (t in THR) {
  cat(sprintf("\n[%s]\n", t))
  if (paste0("SynG_snps_all_", t) %in% names(per))
    cat(sprintf("  significant SNPs, all:     Syn-PALM %s  Obs-GWAS %s  ratio %.2f\n",
                format(s(paste0("SynG_snps_all_", t)), big.mark = ","), format(s(paste0("ObsG_snps_all_", t)), big.mark = ","),
                s(paste0("SynG_snps_all_", t)) / s(paste0("ObsG_snps_all_", t))))
  cat(sprintf("  significant SNPs, no MHC:  Syn-PALM %s  Obs-GWAS %s  ratio %.2f\n",
              format(s(paste0("SynG_snps_noMHC_", t)), big.mark = ","), format(s(paste0("ObsG_snps_noMHC_", t)), big.mark = ","),
              s(paste0("SynG_snps_noMHC_", t)) / s(paste0("ObsG_snps_noMHC_", t))))
  for (m in c("SynG", "ObsG")) {
    z <- lead[[paste(m, t)]]
    if (!nrow(z)) next
    cat(sprintf("  %s loci %s (cis %s, trans %s, gene unmapped %s); proteins with >=1 locus %d\n",
                m, format(nrow(z), big.mark = ","), format(sum(z$cis_trans == "cis", na.rm = TRUE), big.mark = ","),
                format(sum(z$cis_trans == "trans", na.rm = TRUE), big.mark = ","),
                sum(is.na(z$cis_trans)), uniqueN(z$protein)))
  }
}

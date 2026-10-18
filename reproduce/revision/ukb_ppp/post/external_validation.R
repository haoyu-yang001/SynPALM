## External and internal comparison of the primary associations (lead SNPs of
## +-1 Mb clumps at 5e-8, from post_protein.R) of Syn-PALM and Obs-GWAS
## (Results "external validation" paragraph; Reviewer 3 comment 3).
## Per protein and method:
##   catalog : share of loci with a GWAS Catalog variant reported for the protein
##             within +-1 Mb of the lead SNP (GRCh38). A Catalog record belongs
##             to the protein if the gene symbol appears as a word in traitName
##             or the full protein name in efoTraits (the submitted rule, with
##             word boundaries so that e.g. IL1 does not match IL18).
##   sun     : share of loci whose lead SNP has P < 5e-8 in the combined UKB-PPP
##             GWAS of Sun et al. (2023; N = 51,507, which includes the
##             individuals analysed here, so this is a consistency check, not an
##             independent replication); matched by rsid.
## For Syn-PALM loci without an Obs-GWAS locus within +-1 Mb ("Syn-PALM-only"),
## the Obs-GWAS P value, effect estimate and SE at the lead SNP.
##   Rscript external_validation.R <task> <n_task>     per-protein rows
##   Rscript external_validation.R summarize            pooled summary
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
OUT <- "results/table_results_cf_noPRS_noPC/external/"; dir.create(OUT, showWarnings = FALSE)
args <- commandArgs(TRUE)

if (args[1] != "summarize") {
  task <- as.integer(args[1]); ntask <- as.integer(args[2])
  cat_dt <- fread("GWAS_catalog_analysis/files/GWAS_catalog/merged_GWAS_catalog.tsv", sep = "\t", quote = "")
  loc <- sub(";.*", "", cat_dt$locations)
  ok  <- grepl("^[0-9XY]+:[0-9]+", loc)
  cat_dt[, `:=`(cchr = ifelse(ok, sub("^([0-9XY]+):.*", "\\1", loc), NA_character_),
                cpos = ifelse(ok, suppressWarnings(as.integer(sub("^[0-9XY]+:([0-9]+).*", "\\1", loc))), NA_integer_))]
  pn <- fread("results/table_results/protein_names.csv", header = TRUE); setnames(pn, 1:2, c("id", "full"))
  pn[, key := gsub("-", "_", tolower(sub(":.*", "", id)))]
  sun_dirs <- list.dirs("<UKB_DIR>/gwas_data/olink_mapped_results/combined_summary_results", recursive = FALSE)
  sun_key  <- gsub("-", "_", tolower(sub("_.*", "", basename(sun_dirs))))
  files <- sort(list.files("results/table_results_cf_noPRS_noPC", "_loci\\.rds$", full.names = TRUE))
  files <- files[seq(task, length(files), by = ntask)]
  if (nzchar(Sys.getenv("EXT_PROTEINS")))
    files <- files[sub("_loci\\.rds$", "", basename(files)) %in% strsplit(Sys.getenv("EXT_PROTEINS"), ",")[[1]]]
  for (f in files) {
    L <- readRDS(f); p <- L$protein
    o <- paste0(OUT, p, ".rds"); if (file.exists(o)) next
    sym <- toupper(gsub("_", "-", p)); full <- pn$full[match(p, pn$key)]
    rx <- paste0("(^|[^A-Za-z0-9])", gsub("-", "[-_ ]?", sym), "($|[^A-Za-z0-9])")
    cm <- cat_dt[grepl(rx, traitName, ignore.case = TRUE) |
                 (!is.na(full) & grepl(tolower(full), tolower(efoTraits), fixed = TRUE))]
    cm <- cm[!is.na(cpos)]
    loci <- lapply(c(SynG = "SynG", ObsG = "ObsG"), function(m) L$by_threshold[["5e-08"]][[m]]$loci)
    ## Sun et al. statistics, only for the chromosomes that carry a lead SNP
    sd <- sun_dirs[match(p, sun_key)]
    chrs <- unique(as.character(unlist(lapply(loci, function(lo) lo$chr))))
    sf <- if (!is.na(sd)) list.files(sd, "_chr([0-9]+)_summary_results\\.rds$", full.names = TRUE) else character(0)
    sf <- sf[sub(".*_chr([0-9]+)_summary_results\\.rds$", "\\1", sf) %in% chrs]
    sun <- if (!is.na(sd)) rbindlist(lapply(sf, function(s) readRDS(s)[, .(rsid, LOG10P)])) else NULL
    if (!is.null(sun) && !nrow(sun)) sun <- data.table(rsid = character(0), LOG10P = numeric(0))
    rows <- list()
    for (m in names(loci)) {
      lo <- loci[[m]]; if (!nrow(lo)) next
      other <- loci[[setdiff(names(loci), m)]]
      near <- function(chr, pos, d) if (is.null(d) || !nrow(d)) FALSE else any(as.character(d$chr) == chr & abs(d$POS38 - pos) <= 1e6)
      lo[, `:=`(in_catalog = mapply(function(c, x) near(c, x, data.table(chr = cm$cchr, POS38 = cm$cpos)), as.character(chr), POS38),
                sun_sig = if (!is.null(sun)) sun$LOG10P[match(rsid, sun$rsid)] >= -log10(5e-8) else NA,
                only = !mapply(function(c, x) near(c, x, other), as.character(chr), POS38))]
      rows[[m]] <- lo[, .(protein = p, method = m, rsid, chr, POS38, only, in_catalog, sun_sig,
                          SynSurrG_nlp = SynSurrG_negative_log10_pvalues, ObsG_nlp = ObsG_negative_log10_pvalues,
                          SynSurrG_beta = SynSurrG_hat_beta, ObsG_beta = ObsG_hat_beta,
                          SynSurrG_se = sqrt(SynSurrG_var_hat_beta), ObsG_se = sqrt(ObsG_var_hat_beta),
                          catalog_records = nrow(cm), rho_oof = L$rho_oof)]
    }
    saveRDS(rbindlist(rows), o)
  }
  quit(save = "no")
}

## ---- summary ---------------------------------------------------------------
x <- rbindlist(lapply(list.files(OUT, "\\.rds$", full.names = TRUE), readRDS), fill = TRUE)
q3 <- function(v) sprintf("mean %.1f%% (Q1 %.1f%%, Q3 %.1f%%)", 100 * mean(v), 100 * quantile(v, .25), 100 * quantile(v, .75))
pp <- x[catalog_records > 0, .(catalog = mean(in_catalog), n = .N), by = .(protein, method)]
cat("proteins with Catalog records and >= 1 locus:", uniqueN(pp$protein), "\n")
for (m in c("SynG", "ObsG")) cat(sprintf("  %s loci reported in the Catalog (per protein): %s; pooled %.1f%%\n", m,
                                          q3(pp[method == m, catalog]), 100 * x[catalog_records > 0 & method == m, mean(in_catalog)]))
ps <- x[!is.na(sun_sig), .(sun = mean(sun_sig)), by = .(protein, method)]
for (m in c("SynG", "ObsG")) cat(sprintf("  %s loci with P < 5e-8 in Sun et al. (per protein): %s; pooled %.1f%%\n", m,
                                          q3(ps[method == m, sun]), 100 * x[!is.na(sun_sig) & method == m, mean(sun_sig)]))
so <- x[method == "SynG" & only == TRUE]
cat(sprintf("\nSyn-PALM-only loci: %d of %d Syn-PALM loci (%.1f%%)\n", nrow(so), x[method == "SynG", .N], 100 * nrow(so) / x[method == "SynG", .N]))
cat(sprintf("  Obs-GWAS P at the lead SNP: < 0.05 %.1f%%, < 1e-3 %.1f%%, < 1e-5 %.1f%%; median Obs-GWAS -log10 P %.2f\n",
            100 * so[, mean(ObsG_nlp >= -log10(0.05))], 100 * so[, mean(ObsG_nlp >= 3)], 100 * so[, mean(ObsG_nlp >= 5)], so[, median(ObsG_nlp)]))
cat(sprintf("  same sign of the effect estimate %.1f%%; median SE ratio Syn-PALM / Obs-GWAS %.3f; median beta ratio %.3f\n",
            100 * so[, mean(sign(SynSurrG_beta) == sign(ObsG_beta))], so[, median(SynSurrG_se / ObsG_se)], so[, median(SynSurrG_beta / ObsG_beta)]))
cat(sprintf("  in the Catalog %.1f%%; P < 5e-8 in Sun et al. %.1f%%\n", 100 * so[catalog_records > 0, mean(in_catalog)], 100 * so[, mean(sun_sig, na.rm = TRUE)]))
fwrite(x, "results/table_results_cf_noPRS_noPC/external_validation_loci.csv")

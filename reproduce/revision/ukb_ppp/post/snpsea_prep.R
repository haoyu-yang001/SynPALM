## Input SNP sets for SNPsea (Fig 4F, Supplementary Figure obssea; Methods "Expression
## enrichment analysis"): the distinct lead variants of the primary associations
## (P < 5e-8, all proteins) of Syn-PALM, of Obs-GWAS, and of the Syn-PALM associations
## without an Obs-GWAS association within +-1 Mb ("Syn-PALM only"), with GRCh37
## positions (SNPsea data are on hg19) and the smallest P over proteins.
##   Rscript code/post_cf3/snpsea_prep.R ; then code/sh_run/run_noPRS_noPC/run_snpsea.sh
## Writes results/table_results_cf_noPRS_noPC/snpsea/{SynG,ObsG,SynOnly}.gwas
suppressPackageStartupMessages({ library(data.table) })
setwd("<PROJECT_DIR>")
D <- "results/table_results_cf_noPRS_noPC/"; OUT <- paste0(D, "snpsea/"); dir.create(OUT, showWarnings = FALSE)
ann <- as.data.table(readRDS("data/snp_annot_v2geno.rds"))[, .(rsid, POS37)]
mk <- function(z, nlp, f) {
  z <- z[, .(nlp = max(get(nlp))), by = .(rsid, chr)][ann, on = "rsid", nomatch = 0]
  out <- z[, .(CHR = paste0("chr", chr), POS = POS37, SNP = rsid, P = signif(10^(-pmin(nlp, 300)), 3))][order(CHR, POS)]
  fwrite(out, paste0(OUT, f, ".gwas"), sep = "\t"); cat(f, nrow(out), "distinct lead variants\n")
}
mk(fread(paste0(D, "loci_SynG_5e-08.csv")), "SynSurrG_negative_log10_pvalues", "SynG")
mk(fread(paste0(D, "loci_ObsG_5e-08.csv")), "ObsG_negative_log10_pvalues", "ObsG")
ev <- fread(paste0(D, "external_validation_loci.csv"))[method == "SynG" & only == TRUE]
mk(ev[, .(rsid, chr, SynSurrG_negative_log10_pvalues = SynSurrG_nlp)], "SynSurrG_negative_log10_pvalues", "SynOnly")

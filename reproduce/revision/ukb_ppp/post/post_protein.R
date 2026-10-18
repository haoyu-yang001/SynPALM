## Post-processing of one protein's whole-genome scan (no-PRS / no-PC run):
## attach positions, count genome-wide significant SNPs and clump them into
## primary associations exactly as in the submitted analysis (+-1 Mb windows
## around significant SNPs merged, MHC chr6:25.5-34.0 Mb removed, sentinel =
## lowest p), at 5e-8, 1e-8, 1e-9 and 1e-10, for Syn-PALM and Obs-GWAS.
## Autosomes only (the scan also covers chr X, which the paper does not analyse).
## Also writes the protein's genome-wide summary statistics in the web portal's
## format (same file name and columns as hyyang_hufeng/SynGWAS_summary_results/:
## autosomal SNPs with a GRCh38 position, in .bim order) to
## results/web_SynGWAS_summary_results_cf_noPRS_noPC/. Nothing is written to
## the portal's own directory.
##   Rscript post_protein.R <line of order file>
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges) })
setwd("<PROJECT_DIR>")
TAG  <- Sys.getenv("SYNPALM_S1_TAG", "_cf_noPRS_noPC")
line <- as.integer(commandArgs(TRUE)[1])
hdr  <- names(fread("data/df_allUKB_eur1.tsv", nrows = 0)); protein_names <- hdr[49:2968]
p    <- protein_names[as.integer(readLines("code/sh_run/order_cor_vec7_new.txt")[line])]
d    <- paste0("results/allukb_protein_pQTL_european_discovery/", p, "/")
scan <- readRDS(paste0(d, p, TAG, "_scan_allSNPs.rds"))
s1   <- readRDS(paste0(d, p, "_synthetic_phenotype", TAG, ".rds"))
ann  <- readRDS("data/snp_annot_v2geno.rds")
stopifnot(identical(rownames(scan), ann$rsid))
X <- cbind(ann[, .(rsid, chr, POS37, POS38, A1, A2, REF, ALT)], as.data.table(scan))
## p-values that underflow to 0 in -log10 form: recompute from beta / se
for (m in c("SynSurrG", "ObsG", "SynSurr", "Obs")) {
  lp <- paste0(m, "_negative_log10_pvalues")
  z2 <- X[[paste0(m, "_hat_beta")]]^2 / X[[paste0(m, "_var_hat_beta")]]
  bad <- which(!is.finite(X[[lp]]))
  if (length(bad)) set(X, bad, lp, -pchisq(z2[bad], 1, lower.tail = FALSE, log.p = TRUE) / log(10))
}
## (no _withPOS copy: it would duplicate the 57 MB scan for every protein;
##  positions come from data/snp_annot_v2geno.rds whenever needed)

auto <- X$chr %in% as.character(1:22)
keep <- auto & !is.na(X$POS38) & !(X$chr == 6 & X$POS38 >= 25.5e6 & X$POS38 <= 34e6)

## web portal summary statistics
WEB <- "results/web_SynGWAS_summary_results_cf_noPRS_noPC/"
dir.create(WEB, showWarnings = FALSE)
web <- X[auto & !is.na(POS38), .(rsid, chr = as.integer(chr), POS38, REF, ALT,
          SynSurrG_negative_log10_pvalues, SynSurrG_hat_beta, SynSurrG_var_hat_beta,
          ObsG_negative_log10_pvalues, ObsG_hat_beta, ObsG_var_hat_beta)]
stopifnot(!is.unsorted(web$chr))   # .bim (GRCh37) order, as in the portal files
saveRDS(web, paste0(WEB, p, "_negative_log10_pvals_v2geno_crossfit.rds"))
clump <- function(m, thr) {
  lp <- paste0(m, "_negative_log10_pvalues")
  s  <- X[keep & X[[lp]] >= -log10(thr)]
  if (!nrow(s)) return(list(n_snps = 0L, loci = s[0]))
  gr <- GRanges(as.character(s$chr), IRanges(s$POS38 - 1e6, s$POS38 + 1e6))
  rg <- reduce(gr); h <- findOverlaps(gr, rg)
  s[, region := subjectHits(h)]
  lead <- s[s[, .I[which.max(get(lp))], by = region]$V1]
  list(n_snps = nrow(s), loci = lead)
}
out <- list(protein = p, rho_oof = s1$rho_oof, n_obs = s1$n_obs, n_total = attr(scan, "n_total"),
            vc = attr(scan, "vc"), tag = TAG,
            ## all significant autosomal SNPs, MHC and unmapped SNPs included
            n_snps_all = sapply(c(SynG = "SynSurrG", ObsG = "ObsG", SynSurr = "SynSurr", Obs = "Obs"),
                                function(m) sapply(c(5e-8, 1e-8, 1e-9, 1e-10), function(t)
                                  sum(auto & X[[paste0(m, "_negative_log10_pvalues")]] >= -log10(t), na.rm = TRUE))))
for (thr in c(5e-8, 1e-8, 1e-9, 1e-10)) {
  key <- format(thr, scientific = TRUE)
  out$by_threshold[[key]] <- list(SynG = clump("SynSurrG", thr), ObsG = clump("ObsG", thr))
}
dir.create("results/table_results_cf_noPRS_noPC", showWarnings = FALSE)
saveRDS(out, paste0("results/table_results_cf_noPRS_noPC/", p, "_loci.rds"))
b <- out$by_threshold[["5e-08"]]
cat(sprintf("%s rho=%.3f  5e-8: SNPs Syn %d Obs %d | loci Syn %d Obs %d\n", p, s1$rho_oof,
            b$SynG$n_snps, b$ObsG$n_snps, nrow(b$SynG$loci), nrow(b$ObsG$loci)))

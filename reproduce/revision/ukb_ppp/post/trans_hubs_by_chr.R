## Reply to Reviewer 1, Comment 9 (revision, 2026-10-19): do the trans-pQTL signals of
## the UKB-PPP rerun converge on chromosomes 1 and 6, and if so where? Uses the trans
## lead variants and loci of trans_hubs.R (Syn-PALM, P < 1e-10, 250 kb loci, MHC
## chr6:25.5-34.0 Mb excluded).
##   per chromosome: trans signals, distinct proteins, loci with >= 20 / >= 50 proteins,
##     and the share of signals against the share of the tested genotyped autosomal
##     variants (outside the MHC) on the chromosome (obs / exp)
##   chromosomes 1 and 6: every locus with >= 20 proteins (cytoband, lead variant,
##     nearest gene, proteins, per-locus and max-based permutation P of trans_hubs.R)
##   Rscript code/post_cf3/trans_hubs_by_chr.R
## Writes results/table_results_cf_noPRS_noPC/trans_hubs/trans_by_chr_SynG.csv and
## trans_hubs_chr1_chr6_SynG.csv
suppressPackageStartupMessages(library(data.table))
setwd("<PROJECT_DIR>")
D <- "results/table_results_cf_noPRS_noPC/trans_hubs/"
sig <- fread(paste0(D, "trans_signals_SynG.csv"))
hub <- fread(paste0(D, "trans_hubs_SynG.csv"))
ann <- as.data.table(readRDS("data/snp_annot_v2geno.rds"))
ann <- ann[chr %in% as.character(1:22) & !is.na(POS38)]
ann <- ann[!(chr == "6" & POS38 >= 25.5e6 & POS38 <= 34e6)]
nv <- ann[, .(n_variants = .N), by = .(chr = as.integer(chr))]
lc <- sig[, .(n_prot = uniqueN(protein)), by = .(chr, locus_250kb)]
bc <- merge(sig[, .(signals = .N, proteins = uniqueN(protein), loci = uniqueN(locus_250kb)), by = chr],
            lc[, .(loci_ge20 = sum(n_prot >= 20), loci_ge50 = sum(n_prot >= 50), max_proteins = max(n_prot)), by = chr], by = "chr")
bc <- merge(bc, nv, by = "chr")
bc[, `:=`(share_signals = signals / sum(signals), share_variants = n_variants / sum(n_variants))]
bc[, obs_exp := share_signals / share_variants]
bc[, rank_obs_exp := frank(-obs_exp, ties.method = "min")]
setorder(bc, chr)
fwrite(bc, paste0(D, "trans_by_chr_SynG.csv"))
print(bc[, .(chr, signals, proteins, loci, loci_ge20, loci_ge50, max_proteins, share_signals = round(100 * share_signals, 1),
             share_variants = round(100 * share_variants, 1), obs_exp = round(obs_exp, 2), rank_obs_exp)])
cat(sprintf("total trans signals %d; chr1+6 share %.1f%% vs variants %.1f%%\n", sum(bc$signals),
            100 * bc[chr %in% c(1, 6), sum(share_signals)], 100 * bc[chr %in% c(1, 6), sum(share_variants)]))
h16 <- hub[chr %in% c(1, 6) & n_proteins >= 20,
           .(chr, interval_Mb = sprintf("%.2f-%.2f", start_Mb, end_Mb), cytoband, lead_rsid, nearest_gene, n_proteins,
             perm_P_locus, perm_P_max, median_abs_r, sign_concord_at_lead, ObsG_n_proteins)][order(chr, -n_proteins)]
fwrite(h16, paste0(D, "trans_hubs_chr1_chr6_SynG.csv"))
print(h16)
cat(sprintf("loci with >= 20 proteins: genome %d; chr1 %d, chr6 %d; genome-wide rank of the largest chr1 / chr6 locus: %d / %d\n",
            nrow(hub[n_proteins >= 20]), nrow(h16[chr == 1]), nrow(h16[chr == 6]),
            match(TRUE, hub[order(-n_proteins)]$chr == 1), match(TRUE, hub[order(-n_proteins)]$chr == 6)))

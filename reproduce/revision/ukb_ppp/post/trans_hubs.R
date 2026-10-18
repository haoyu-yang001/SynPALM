## Candidate pleiotropic trans loci ("trans hubs") of the no-PRS / no-PC rerun
## (UKB-PPP), from the P < 1e-10 locus lead variants written by aggregate_all.R;
## Syn-PALM is the main analysis, Obs-GWAS is run identically for comparison.
## Method (as in the manuscript Methods): autosomes only, HLA region
## chr6:25.5-34.0 Mb excluded; a lead variant is trans if it lies on another
## chromosome than, or > 1 Mb outside the body (start_position-end_position of
## final_protein_map_38.rds) of, the protein's coding gene (proteins without gene
## coordinates are dropped); trans lead variants of all proteins are merged by
## greedy positional clustering (gap 250 kb; 500 kb and 1 Mb as sensitivity) and
## the distinct proteins per locus are counted. Permutation (N_PERM): each
## protein keeps its number of trans signals, which are re-drawn (without
## replacement) from the fixed set of distinct observed trans lead positions;
## distinct proteins per locus are recounted and the maximum recorded
## (P for the observed maximum and, max-based, for every locus; a per-locus P is
## also given). A second, multiplicity-preserving null (protein labels shuffled
## across all trans signals) is reported for reference. For loci with >= K
## proteins: most frequent lead rsid, nearest gene (TxDb hg38 knownGene, as
## get_nearest_gene_from_pos in code/gwas_plot_function.R), cytoband (UCSC hg38),
## 4 example proteins with largest |beta| of their own lead variant, median
## |Pearson r| of measured levels among the implicated proteins (pairwise
## complete obs., data/df_allUKB_eur1.tsv), and sign concordance (share of
## proteins with the majority sign) of each protein's effect at the locus's most
## frequent lead rsid (same allele for all proteins; from the per-protein
## summary statistics in results/web_SynGWAS_summary_results_cf_noPRS_noPC/),
## plus, unaligned, of each protein's own lead variant.
##   sbatch -p test -c 8 --mem=64G -t 0-2:00 --wrap \
##     "module load R/4.4.1-fasrc01; Rscript code/post_cf3/trans_hubs.R"
## Writes to results/table_results_cf_noPRS_noPC/trans_hubs/:
##   trans_hubs_SynG.csv, trans_hubs_ObsG.csv  one row per trans locus with >= K
##       proteins (250 kb gap), ranked by n_proteins; the Syn-PALM file also has
##       the overlapping Obs-GWAS locus (n proteins, lead rsid, shared proteins)
##   trans_signals_SynG.csv, trans_signals_ObsG.csv  every trans lead variant
##       with its locus (250 kb / 500 kb / 1 Mb gap)
##   trans_hubs_summary.txt  numbers for the manuscript paragraph
##   cytoBand_hg38.txt.gz    UCSC cytoband table (downloaded once, then cached)
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(parallel) })
setwd("<PROJECT_DIR>")
D    <- "results/table_results_cf_noPRS_noPC/"
OUT  <- paste0(D, "trans_hubs/"); dir.create(OUT, showWarnings = FALSE)
WEB  <- "results/web_SynGWAS_summary_results_cf_noPRS_noPC/"
THR  <- "1e-10"
CIS  <- 1e6                       # cis window around the gene body
GAPS <- c(2.5e5, 5e5, 1e6)        # first = main
HLA  <- c(25.5e6, 34e6)           # chr6, excluded
K    <- 20                        # loci with >= K proteins are tabulated
N_TOP <- 8; N_EX <- 4; N_PERM <- 200; SEED <- 20261013
## named loci quoted in the manuscript, looked up by rsid (reported if present)
LOI  <- c("rs61747728 (1q25.2, NPHS2)" = "rs61747728", "rs7775698 (6q23.3, HBS1L-MYB)" = "rs7775698")
NCORE <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
BCOL <- c(SynG = "SynSurrG_hat_beta", ObsG = "ObsG_hat_beta")
PCOL <- c(SynG = "SynSurrG_negative_log10_pvalues", ObsG = "ObsG_negative_log10_pvalues")
LAB  <- c(SynG = "Syn-PALM", ObsG = "Obs-GWAS")
f <- function(x) format(x, big.mark = ",")

## ---- gene coordinates ---------------------------------------------------
gmap <- as.data.table(readRDS("<LAB_DIR>/protein_map/final_protein_map_38.rds"))
gmap <- unique(gmap[, .(protein = gsub("-", "_", id), gene = gene_upper, gene_chr = as.character(chr),
                        gene_start = start_position, gene_end = end_position)], by = "protein")
gmap <- gmap[!is.na(gene_chr) & !is.na(gene_start) & !is.na(gene_end)]
analysed <- fread(paste0(D, "summary_per_protein.csv"), select = "protein")$protein

## ---- trans lead variants ------------------------------------------------
read_trans <- function(m) {
  z <- fread(paste0(D, "loci_", m, "_", THR, ".csv"))
  z <- z[, .(protein, rsid, chr = as.integer(chr), POS38, beta = get(BCOL[m]), lp = get(PCOL[m]))]
  z <- z[chr %in% 1:22 & !is.na(POS38) & !(chr == 6 & POS38 >= HLA[1] & POS38 <= HLA[2])]
  all_n <- c(signals = nrow(z), proteins = uniqueN(z$protein))
  z <- merge(z, gmap, by = "protein")
  z[, dist := fifelse(as.character(chr) != gene_chr, Inf, pmax(0, gene_start - POS38, POS38 - gene_end))]
  attr(z, "n_all") <- all_n
  z
}

## greedy positional clustering of the distinct positions, per chromosome
cluster <- function(chr, pos, gap) {
  q <- data.table(chr, pos)
  u <- unique(q)[order(chr, pos)]
  u[, cl := cumsum(c(TRUE, diff(pos) > gap)), by = chr]
  u[, locus := sprintf("chr%d_%03d", chr, cl)]
  u[q, on = .(chr, pos), locus]
}

count_loci <- function(locus, protein) {
  x <- unique(data.table(locus, protein))[, .N, by = locus]
  setNames(x$N, x$locus)
}

## null: proteins keep their number of trans signals, positions re-drawn from
## the set of distinct observed positions; returns n_perm x n_loci counts
perm_counts <- function(t, loci_names, n_perm) {
  pool <- unique(t[, .(chr, POS38, locus)])
  k <- t[, .(k = uniqueN(paste(chr, POS38))), by = protein]$k
  lab <- rep(seq_along(k), k)
  t(vapply(seq_len(n_perm), function(i) {
    idx <- unlist(lapply(k, function(kk) sample.int(nrow(pool), kk)), use.names = FALSE)
    cnt <- count_loci(pool$locus[idx], lab)
    out <- cnt[loci_names]; out[is.na(out)] <- 0L; out
  }, numeric(length(loci_names))))
}
## reference null: protein labels shuffled across all trans signals
perm_shuffle_max <- function(t, n_perm)
  vapply(seq_len(n_perm), function(i) max(count_loci(t$locus, sample(t$protein))), 0)

## ---- annotation helpers ---------------------------------------------------
cyto_file <- paste0(OUT, "cytoBand_hg38.txt.gz")
if (!file.exists(cyto_file))
  try(download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
                    cyto_file, quiet = TRUE), silent = TRUE)
cyto <- if (file.exists(cyto_file)) fread(cyto_file, col.names = c("chr", "start", "end", "band", "stain")) else NULL
band_of <- function(cc, s0, e0) {
  if (is.null(cyto)) return(NA_character_)
  b <- cyto[chr == paste0("chr", cc) & end >= s0 & start < e0, band]   # UCSC 0-based half-open
  if (!length(b)) return(NA_character_)
  if (length(b) == 1) paste0(cc, b[1]) else paste0(cc, b[1], "-", sub("^[pq]", "", b[length(b)]))
}
nearest_gene <- function(chr, pos) {   # as get_nearest_gene_from_pos(build = "hg38"), named genes only
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  gl <- suppressMessages(GenomicFeatures::genes(txdb, single.strand.genes.only = FALSE))
  g1 <- unlist(gl); g1 <- g1[as.character(seqnames(g1)) %in% paste0("chr", c(1:22, "X", "Y"))]
  gg <- unlist(range(split(g1, names(g1))))   # genes with alt-contig copies (e.g. ABO) are kept
  gg$symbol <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(gg),
                                                             "SYMBOL", "ENTREZID")))
  ## unnamed / non-coding-style symbols (LOC, LINC, MIR, SNOR, antisense) are skipped
  gg <- gg[!is.na(gg$symbol) & !grepl("^(LOC|LINC|MIR|SNOR)|-AS[0-9]*$|-DT$", gg$symbol)]
  q <- GRanges(paste0("chr", chr), IRanges(pos, pos))
  h <- nearest(q, gg, ignore.strand = TRUE, select = "all")
  d <- distance(q[queryHits(h)], gg[subjectHits(h)], ignore.strand = TRUE)
  x <- data.table(i = queryHits(h), g = gg$symbol[subjectHits(h)], d = d)
  x <- x[, .(nearest_gene = paste(sort(unique(g)), collapse = "/"), dist_to_gene_kb = round(d[1] / 1e3, 1)), by = i]
  x[data.table(i = seq_along(q)), on = "i"][, .(nearest_gene, dist_to_gene_kb)]
}

## ---- run one method ---------------------------------------------------------
run <- function(m) {
  set.seed(SEED)
  z <- read_trans(m)
  t <- z[dist > CIS]
  for (g in GAPS) t[, paste0("locus_", g / 1e3, "kb") := cluster(chr, POS38, g)]
  t[, locus := get(paste0("locus_", GAPS[1] / 1e3, "kb"))]
  ## locus summaries for every gap, with permutation
  by_gap <- lapply(GAPS, function(g) {
    tt <- copy(t)[, locus := get(paste0("locus_", g / 1e3, "kb"))]
    L <- tt[, .(chr = chr[1], start = min(POS38), end = max(POS38), n_signals = .N,
                n_proteins = uniqueN(protein), n_variants = uniqueN(rsid)), by = locus]
    setorder(L, -n_proteins, chr, start)
    P <- perm_counts(tt, L$locus, N_PERM)
    nmax <- apply(P, 1, max)
    L[, perm_P_max := (1 + vapply(n_proteins, function(n) sum(nmax >= n), 0)) / (N_PERM + 1)]
    L[, perm_P_locus := (1 + colSums(sweep(P, 2, n_proteins, ">="))) / (N_PERM + 1)]
    L[, perm_mean_locus := colMeans(P)]
    list(L = L, nmax = nmax, shuf = perm_shuffle_max(tt, N_PERM))
  })
  names(by_gap) <- paste0(GAPS / 1e3, "kb")
  L <- by_gap[[1]]$L
  ## lead rsid = most frequent lead rsid in the locus (ties: strongest signal)
  lead <- t[, .(n_prot_at_lead = uniqueN(protein), maxlp = max(lp), lead_pos = POS38[1]), by = .(locus, rsid)]
  setorder(lead, locus, -n_prot_at_lead, -maxlp)
  lead <- lead[, .SD[1], by = locus][, .(locus, lead_rsid = rsid, lead_pos, n_prot_at_lead)]
  L <- merge(L, lead, by = "locus", sort = FALSE)
  ex <- t[, .(ab = max(abs(beta))), by = .(locus, protein)][order(locus, -ab)]
  ex <- ex[, .(example_proteins = paste(toupper(head(protein, N_EX)), collapse = ", "),
               proteins = paste(protein, collapse = ";")), by = locus]
  L <- merge(L, ex, by = "locus", sort = FALSE)
  setorder(L, -n_proteins, chr, start)
  list(m = m, z = z, t = t, L = L, by_gap = by_gap, n_all = attr(z, "n_all"))
}
msg <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)
msg("loci and permutations")
R <- lapply(c(SynG = "SynG", ObsG = "ObsG"), run)
msg("annotation")

## ---- annotate loci with >= K proteins -----------------------------------
for (m in names(R)) {
  H <- R[[m]]$L[n_proteins >= K]
  H[, cytoband := mapply(band_of, chr, start, end)]
  H[, cytoband_lead := mapply(band_of, chr, lead_pos, lead_pos)]
  H <- cbind(H, nearest_gene(H$chr, H$lead_pos))
  R[[m]]$H <- H
}

## effects of every implicated protein at its locus's lead rsid (same allele)
need <- unique(rbindlist(lapply(R, function(r)
  r$H[, .(protein = strsplit(proteins, ";")[[1]]), by = .(lead_rsid)])))
msg("effects at lead rsids")
eff <- rbindlist(mclapply(unique(need$protein), function(p) {
  fn <- paste0(WEB, p, "_negative_log10_pvals_v2geno_crossfit.rds")
  if (!file.exists(fn)) return(NULL)            # protein then counts as missing in n_prot_sign
  w <- readRDS(fn)
  w[rsid %in% need[protein == p]$lead_rsid,
    .(protein = p, lead_rsid = rsid, SynG = SynSurrG_hat_beta, ObsG = ObsG_hat_beta)]
}, mc.cores = NCORE, mc.preschedule = FALSE))
no_sumstats <- setdiff(unique(need$protein), eff$protein)
conc <- function(b) { b <- b[is.finite(b) & b != 0]; if (!length(b)) NA_real_ else max(mean(b > 0), mean(b < 0)) }

## measured protein levels
allp <- unique(unlist(lapply(R, function(r) strsplit(r$H$proteins, ";"))))
msg("measured levels")
hdr  <- names(fread("data/df_allUKB_eur1.tsv", nrows = 0))
lev  <- fread("data/df_allUKB_eur1.tsv", select = intersect(allp, hdr), sep = "\t")
lev  <- lev[rowSums(!is.na(lev)) > 0]          # participants with any measured protein
## = cor(X, use = "pairwise.complete.obs"), via cross-products (fast for many columns)
pcor <- function(X) {
  X <- sweep(X, 2, colMeans(X, na.rm = TRUE)); M <- 1 * !is.na(X); X[is.na(X)] <- 0
  n <- crossprod(M); sx <- crossprod(X, M); sxx <- crossprod(X^2, M); v <- n * sxx - sx^2
  (n * crossprod(X) - sx * t(sx)) / sqrt(v * t(v))
}
cor_med <- function(ps) {
  ps <- intersect(ps, names(lev)); if (length(ps) < 2) return(c(NA, length(ps)))
  C <- pcor(as.matrix(lev[, ..ps]))
  c(median(abs(C[lower.tri(C)]), na.rm = TRUE), length(ps))
}

msg("correlations, sign concordance")
for (m in names(R)) {
  H <- R[[m]]$H; t <- R[[m]]$t
  st <- rbindlist(mclapply(seq_len(nrow(H)), function(i) {
    ps <- strsplit(H$proteins[i], ";")[[1]]
    b_lead <- eff[lead_rsid == H$lead_rsid[i] & protein %in% ps][[m]]
    own <- t[locus == H$locus[i]][order(-lp)][, .SD[1], by = protein]$beta
    cm <- cor_med(ps)
    data.table(median_abs_r = cm[1], n_prot_r = cm[2],
               sign_concord_at_lead = conc(b_lead), n_prot_sign = sum(is.finite(b_lead)),
               sign_concord_own_lead_unaligned = conc(own))
  }, mc.cores = NCORE))
  R[[m]]$H <- cbind(H, st)
}

## Obs-GWAS locus overlapping each Syn-PALM locus (+- main gap)
HS <- R$SynG$H; LO <- R$ObsG$L
cmp <- rbindlist(lapply(seq_len(nrow(HS)), function(i) {
  o <- LO[chr == HS$chr[i] & end >= HS$start[i] - GAPS[1] & start <= HS$end[i] + GAPS[1]][order(-n_proteins)]
  if (!nrow(o)) return(data.table(ObsG_locus = NA_character_, ObsG_n_proteins = 0L, ObsG_lead_rsid = NA_character_,
                                  ObsG_n_shared_proteins = 0L))
  sp <- intersect(strsplit(HS$proteins[i], ";")[[1]], unlist(strsplit(o$proteins, ";")))
  data.table(ObsG_locus = paste(o$locus, collapse = ";"), ObsG_n_proteins = o$n_proteins[1],
             ObsG_lead_rsid = o$lead_rsid[1], ObsG_n_shared_proteins = length(sp))
}))
R$SynG$H <- cbind(HS, cmp)

## ---- write --------------------------------------------------------------
cols <- c("locus", "chr", "start", "end", "cytoband", "cytoband_lead", "lead_rsid", "lead_pos", "n_prot_at_lead",
          "nearest_gene", "dist_to_gene_kb", "n_signals", "n_variants", "n_proteins", "perm_P_max", "perm_P_locus",
          "perm_mean_locus", "example_proteins", "median_abs_r", "n_prot_r", "sign_concord_at_lead", "n_prot_sign",
          "sign_concord_own_lead_unaligned")
for (m in names(R)) {
  H <- R[[m]]$H
  H[, `:=`(start_Mb = round(start / 1e6, 3), end_Mb = round(end / 1e6, 3))]
  keep <- c(sub("^start$", "start_Mb", sub("^end$", "end_Mb", cols)), grep("^ObsG_", names(H), value = TRUE), "proteins")
  fwrite(H[, ..keep], paste0(OUT, "trans_hubs_", m, ".csv"))
  fwrite(R[[m]]$t[, c("protein", "gene", "gene_chr", "gene_start", "gene_end", "rsid", "chr", "POS38", "beta", "lp",
                      grep("^locus_", names(R[[m]]$t), value = TRUE)), with = FALSE],
         paste0(OUT, "trans_signals_", m, ".csv"))
}

## ---- summary --------------------------------------------------------------
options(width = 250)
sink(paste0(OUT, "trans_hubs_summary.txt"), split = TRUE)
cat(sprintf("Candidate pleiotropic trans loci, P < %s locus lead variants (aggregate_all.R), autosomes,\n", THR))
cat(sprintf("HLA chr6:%.1f-%.1f Mb excluded; trans = other chromosome or > %g Mb outside the gene body;\n",
            HLA[1] / 1e6, HLA[2] / 1e6, CIS / 1e6))
cat(sprintf("greedy clustering gap %s kb (main); permutation n = %d, seed %d. Run %s\n\n",
            GAPS[1] / 1e3, N_PERM, SEED, format(Sys.time(), "%Y-%m-%d %H:%M")))
cat(sprintf("proteins analysed: %s; with gene coordinates: %s (gene on an autosome: %s, chr X: %s)\n",
            f(length(analysed)), f(sum(analysed %in% gmap$protein)),
            f(sum(analysed %in% gmap[gene_chr %in% 1:22]$protein)), f(sum(analysed %in% gmap[gene_chr == "X"]$protein))))
for (m in names(R)) {
  r <- R[[m]]; L <- r$L; t <- r$t; z <- r$z; bg <- r$by_gap
  cat(sprintf("\n======== %s ========\n", LAB[m]))
  cat(sprintf("lead variants (autosomal, non-HLA): %s for %s proteins; with gene coordinates: %s for %s proteins\n",
              f(r$n_all["signals"]), f(r$n_all["proteins"]), f(nrow(z)), f(uniqueN(z$protein))))
  cat(sprintf("cis %s, trans %s (trans on another chromosome: %s)\n", f(sum(z$dist <= CIS)), f(nrow(t)), f(sum(is.infinite(t$dist)))))
  cat(sprintf("TRANS: %s signals for %s proteins across %s loci (%s distinct lead variants)\n",
              f(nrow(t)), f(uniqueN(t$protein)), f(nrow(L)), f(uniqueN(t$rsid))))
  cat(sprintf("proteins per locus: median %g, IQR %g-%g, max %d; loci with >=5 / >=%d / >=50 / >=100 proteins: %d / %d / %d / %d\n",
              median(L$n_proteins), quantile(L$n_proteins, .25), quantile(L$n_proteins, .75), max(L$n_proteins),
              K, sum(L$n_proteins >= 5), sum(L$n_proteins >= K), sum(L$n_proteins >= 50), sum(L$n_proteins >= 100)))
  cat("\nlocus-width sensitivity and permutation (null max of distinct proteins per locus):\n")
  print(as.data.frame(rbindlist(lapply(names(bg), function(g) { x <- bg[[g]]
    data.table(gap = g, n_loci = nrow(x$L), median = median(x$L$n_proteins), n_ge50 = sum(x$L$n_proteins >= 50),
               max_obs = max(x$L$n_proteins), null_max_mean = mean(x$nmax), null_max_q95 = quantile(x$nmax, .95),
               null_max_max = max(x$nmax), P_obs_max = (1 + sum(x$nmax >= max(x$L$n_proteins))) / (N_PERM + 1),
               loci_P_lt_0.05 = sum(x$L$perm_P_max < 0.05),
               shuffle_null_max_mean = mean(x$shuf), shuffle_P = (1 + sum(x$shuf >= max(x$L$n_proteins))) / (N_PERM + 1))
  }))), row.names = FALSE)
  cat("(perm: each protein keeps its number of trans signals, positions drawn from the set of distinct observed\n",
      " trans lead positions; shuffle: protein labels permuted across signals, positions kept with multiplicity)\n", sep = "")
  H <- r$H
  cat(sprintf("\nloci with >= %d proteins: %d; top %d:\n", K, nrow(H), N_TOP))
  sh <- head(H, N_TOP)[, .(rank = .I, locus, interval_Mb = sprintf("%d:%.2f-%.2f", chr, start / 1e6, end / 1e6),
                           cytoband, nearest_gene, lead_rsid, n_at_lead = n_prot_at_lead, n_proteins,
                           P_max = signif(perm_P_max, 2), med_abs_r = round(median_abs_r, 3),
                           sign_at_lead = round(sign_concord_at_lead, 2), sign_own = round(sign_concord_own_lead_unaligned, 2))]
  print(as.data.frame(sh), row.names = FALSE)
  cat("examples (largest |beta| of own lead variant):\n")
  for (i in seq_len(min(N_TOP, nrow(H)))) cat(sprintf("  %d %s %s: %s\n", i, H$cytoband[i], H$nearest_gene[i], H$example_proteins[i]))
  cat(sprintf("top-%d medians: median |r| %.3f, sign concordance at lead rsid %.2f\n", N_TOP,
              median(head(H$median_abs_r, N_TOP)), median(head(H$sign_concord_at_lead, N_TOP))))
  if (m == "SynG") {
    cat("\nObs-GWAS locus overlapping each top Syn-PALM locus:\n")
    print(as.data.frame(head(H, N_TOP)[, .(locus, nearest_gene, lead_rsid, n_proteins, ObsG_locus, ObsG_lead_rsid,
                                           ObsG_n_proteins, ObsG_n_shared_proteins)]), row.names = FALSE)
    ot <- head(R$ObsG$H, N_TOP)
    cat(sprintf("top-%d by Obs-GWAS lead rsids: %s\n", N_TOP, paste(ot$lead_rsid, collapse = ", ")))
    cat(sprintf("top-%d lead rsids shared with Syn-PALM top-%d: %d\n", N_TOP, N_TOP,
                length(intersect(ot$lead_rsid, head(H$lead_rsid, N_TOP)))))
  }
  cat("\nnamed loci:\n")
  for (nm in names(LOI)) {
    l <- unique(t[rsid == LOI[nm]]$locus)
    if (!length(l)) { cat(sprintf("  %s: not a trans lead variant\n", nm)); next }
    x <- L[locus %in% l]
    for (j in seq_len(nrow(x)))
      cat(sprintf("  %s: locus %s chr%d:%.2f-%.2f Mb, %d proteins (%d with this lead), locus lead %s, perm P_max %.3g, P_locus %.3g\n",
                  nm, x$locus[j], x$chr[j], x$start[j] / 1e6, x$end[j] / 1e6, x$n_proteins[j],
                  uniqueN(t[rsid == LOI[nm] & locus == x$locus[j]]$protein), x$lead_rsid[j], x$perm_P_max[j], x$perm_P_locus[j]))
  }
}
if (length(no_sumstats)) cat(sprintf("\nno portal summary statistics (excluded from sign concordance at lead rsid): %s\n",
                                     paste(no_sumstats, collapse = ", ")))
cat("\nOutputs:", paste0(OUT, c("trans_hubs_SynG.csv", "trans_hubs_ObsG.csv", "trans_signals_SynG.csv",
                                "trans_signals_ObsG.csv")), sep = "\n  ")
sink()

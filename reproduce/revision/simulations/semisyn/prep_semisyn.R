## =====================================================================
## prep_semisyn.R -- one-off design of the semi-synthetic simulation
## (NC R1 revision; reviewers: no simulation trains the prediction model
## inside the loop).
##
## Real UKB surrogates S (data/extracted_pheno_data528_eur.tsv + ICD-10
## level-1 table), real genotypes (ukb_cal_allchr_v2_whiteb_QC), real sparse
## GRM and the UKB-PPP relatedness-cluster folds (data/crossfit_folds_K5.rds).
## This script fixes everything that does NOT change between replicates:
##
##  1. cluster-level re-assignment of the surrogate rows. The rows of the
##     surrogate + ICD tables are moved, jointly, from one relatedness cluster
##     (connected component of the sparse GRM) to another of the same size and
##     sex pattern; singletons and clusters without a match are permuted one by
##     one within sex (seed SEED_PERM).
##     S keeps its real joint distribution, its real within-family correlation
##     and sex, but is independent of the genotypes. Without this every
##     genotyped SNP would carry a small real S-mediated effect on
##     Y* = h(S) + ..., and no SNP would be an exact null.
##  2. the fixed SNP set (autosomes only):
##       causal_Y   80 SNPs, odd chromosomes, effect on Y only
##       causal_YS  80 SNPs, odd chromosomes, effect on Y and on surrogate Z
##       S_null   2000 SNPs, odd chromosomes, >= 1 Mb from every causal SNP,
##                effect on ONE surrogate column only (400 per column; 4 columns
##                enter h, the 5th is Z) -- null for Y, non-null for Yhat
##       null    30000 SNPs, even chromosomes (no planted SNP there)
##     causal SNPs >= 2 Mb apart, S_null >= 250 kb apart, MAF >= 0.05 and call
##     rate >= 0.99 for planted SNPs, MAF >= 0.01 and call rate >= 0.98 for nulls.
##  3. the planted effects (seed SEED_EFF) and the planted genetic scores,
##     sum_j b_j (G_j - 2p_j) for Y* and for each surrogate column (missing
##     genotypes set to 0 as in the scan, centred with the same 2p_j).
##
## Outputs (cluster only, never copied to git):
##   simulation_code/results_rev1/semisyn/design.rds          (SNP table, effects, seeds, h)
##   simulation_code/results_rev1/semisyn/design_indiv.rds    (row permutation, planted scores:
##                                                             individual level -- stays on the cluster)
##
##   Rscript prep_semisyn.R
## =====================================================================
setwd("<PROJECT_DIR>")
suppressMessages({library(data.table); library(Matrix); library(withr)})
OUT <- "simulation_code/results_rev1/semisyn/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
tick <- function(...) cat(sprintf("[%6.1f s] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

SEED_PERM <- 20261020L; SEED_SNP <- 20261021L; SEED_EFF <- 20261022L
geno_prefix <- "genotype_data/ukb_cal_allchr_v2_whiteb_QC"

## surrogate columns: four enter h(S) (BMI, SBP, HbA1c, LDL), Z (platelet count)
## carries the shared polygenic term and the Y+S causal SNPs but is not in h
H_COLS <- c(bmi = "p21001_i0", sbp = "p4080_i0_a0", hba1c = "p30750_i0", ldl = "p30780_i0")
Z_COL  <- c(plt = "p30080_i0")
S_COLS <- c(H_COLS, Z_COL)

FD <- readRDS("data/crossfit_folds_K5.rds")
n  <- length(FD$analysis_eid)
S  <- fread("data/extracted_pheno_data528_eur.tsv", data.table = FALSE, select = c("eid", "p31", unname(S_COLS)))
S  <- S[match(FD$analysis_eid, S$eid), ]
stopifnot(all(as.numeric(S$eid) == as.numeric(FD$analysis_eid)))

### 1. row re-assignment between relatedness clusters of the same size and sex pattern
## clusters = connected components of the sparse GRM (FD$blocks are contiguous
## index ranges that can join unrelated people, so they are not used here)
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
GRM <- sgrm[FD$index_all, FD$index_all]; rm(sgrm); invisible(gc())
sm <- summary(as(GRM, "generalMatrix")); sm <- sm[sm$i != sm$j & sm$x != 0, ]
lab <- seq_len(n)
repeat {                                   # min-label propagation + pointer jumping
  m <- pmin(lab[sm$i], lab[sm$j]); agg <- tapply(m, sm$i, min)
  new <- lab; ii <- as.integer(names(agg)); new[ii] <- pmin(new[ii], agg); new <- new[new]
  if (identical(new, lab)) break
  lab <- new
}
comps <- unname(split(seq_len(n), lab)); csz <- lengths(comps)
tick("relatedness clusters: ", length(comps), ", multi-member = ", sum(csz > 1), ", largest = ", max(csz),
     ", individuals in clusters: ", sum(csz[csz > 1]))
rm(GRM, sm); invisible(gc())
key <- vapply(comps, function(b) paste0(length(b), ":", paste(S$p31[b], collapse = "")), "")
perm <- seq_len(n)
leftover <- integer(0)
with_seed(SEED_PERM, {
  for (k in unique(key)) {
    ii <- which(key == k)
    if (length(ii) > 1 && csz[ii[1]] > 1) {           # whole clusters to clusters
      don <- ii[sample.int(length(ii))]
      for (t in seq_along(ii)) perm[comps[[ii[t]]]] <- comps[[don[t]]]
    } else leftover <- c(leftover, unlist(comps[ii]))  # singletons and unmatched clusters
  }
  for (sx in unique(S$p31[leftover])) {               # individuals, within sex
    jj <- leftover[S$p31[leftover] == sx]
    perm[jj] <- jj[sample.int(length(jj))]
  }
})
stopifnot(identical(sort(perm), seq_len(n)), all(S$p31[perm] == S$p31))
moved <- sum(perm != seq_len(n))
tick("rows moved: ", moved, " of ", n, "; individuals in unmatched clusters (moved one by one): ",
     sum(csz[csz > 1 & !(key %in% key[duplicated(key)])]))

### 2. SNP set ###############################################################
fam <- fread(paste0(geno_prefix, ".fam"), header = FALSE, data.table = FALSE)
keep <- fam[match(as.character(FD$analysis_eid), as.character(fam$V2)), 1:2]
stopifnot(!anyNA(keep$V2))
kf <- paste0(OUT, "tmp_keep.txt"); fwrite(keep, kf, sep = "\t", col.names = FALSE)
frq_pref <- paste0(OUT, "tmp_freq")
if (!file.exists(paste0(frq_pref, ".afreq")))
  system2("plink2", c("--bfile", geno_prefix, "--keep", kf, "--freq",
                                            "--threads", "4", "--memory", "8000", "--out", frq_pref))
frq <- fread(paste0(frq_pref, ".afreq"), data.table = FALSE)
bim <- fread(paste0(geno_prefix, ".bim"), header = FALSE, data.table = FALSE)
stopifnot(identical(frq$ID, bim$V2))
snp <- data.frame(rsid = bim$V2, chr = suppressWarnings(as.integer(bim$V1)), pos = bim$V4,
                  a1 = bim$V5, alt_freq = frq$ALT_FREQS, call = frq$OBS_CT / (2 * n))
snp$maf <- pmin(snp$alt_freq, 1 - snp$alt_freq)
snp <- snp[!is.na(snp$chr) & snp$chr %in% 1:22, ]
tick("autosomal SNPs: ", nrow(snp))

## greedy random pick with a minimum distance to the already-picked SNPs
pick_spaced <- function(pool, m, min_dist, avoid = NULL, avoid_dist = 0) {
  pool <- pool[sample.int(nrow(pool)), ]
  if (!is.null(avoid)) {
    bad <- logical(nrow(pool))
    for (cc in unique(avoid$chr)) {
      i <- which(pool$chr == cc); a <- avoid$pos[avoid$chr == cc]
      if (length(i)) bad[i] <- vapply(pool$pos[i], function(p) any(abs(p - a) < avoid_dist), logical(1))
    }
    pool <- pool[!bad, ]
  }
  took <- integer(0); last <- list()
  for (r in seq_len(nrow(pool))) {
    cc <- as.character(pool$chr[r]); p <- pool$pos[r]
    if (!is.null(last[[cc]]) && any(abs(p - last[[cc]]) < min_dist)) next
    took <- c(took, r); last[[cc]] <- c(last[[cc]], p)
    if (length(took) == m) break
  }
  if (length(took) < m) stop("only ", length(took), " SNPs could be placed (wanted ", m, ")")
  pool[took, ]
}
odd  <- snp$chr %% 2 == 1
planted_ok <- snp$maf >= 0.05 & snp$call >= 0.99
with_seed(SEED_SNP, {
  causal <- pick_spaced(snp[odd & planted_ok, ], 160, min_dist = 2e6)
  snull  <- pick_spaced(snp[odd & planted_ok & !snp$rsid %in% causal$rsid, ], 2000, min_dist = 2.5e5,
                        avoid = causal, avoid_dist = 1e6)
  nullp  <- snp[!odd & snp$maf >= 0.01 & snp$call >= 0.98, ]
  rnull  <- nullp[sort(sample.int(nrow(nullp), 30000)), ]
})
causal$cat <- rep(c("causal_Y", "causal_YS"), each = 80)
causal$q   <- rep(rep(c(1, 2, 4, 8) * 1e-4, each = 20), 2)       # share of Var(Y*)
snull$cat  <- "S_null"
snull$col  <- rep(names(S_COLS), each = 400)
snull$qS   <- rep(c(rep(4e-4, 390), rep(5e-3, 10)), 5)           # share of Var(S0_col)
rnull$cat  <- "null"
tick("picked: ", nrow(causal), " causal, ", nrow(snull), " S-acting nulls, ", nrow(rnull), " random nulls")

### 3. effects ###############################################################
QZ_YS <- 2.5e-3                                # causal_YS: share of Var(Z0) explained
sdS0  <- vapply(S[, S_COLS], sd, 0); names(sdS0) <- names(S_COLS)
with_seed(SEED_EFF, {
  causal$sign  <- sample(c(-1, 1), nrow(causal), replace = TRUE)
  snull$sign   <- sample(c(-1, 1), nrow(snull),  replace = TRUE)
})
v2pq <- function(d) 2 * d$alt_freq * (1 - d$alt_freq)
causal$b_Y   <- causal$sign * sqrt(causal$q / v2pq(causal))          # per ALT allele, Var(Y*) ~ 1
causal$g_Z   <- ifelse(causal$cat == "causal_YS", causal$sign * sdS0[["plt"]] * sqrt(QZ_YS / v2pq(causal)), 0)
snull$g_S    <- snull$sign * sdS0[snull$col] * sqrt(snull$qS / v2pq(snull))

## planted scores. Effects are per BEDMatrix-counted allele (the dosage the
## scan tests); 2pq does not depend on which allele is counted and the sign is
## random, so the variance shares hold. Centring uses the data frequency.
G0 <- BEDMatrix::BEDMatrix(paste0(geno_prefix, ".bed"), simple_names = TRUE)
rows <- match(as.character(FD$analysis_eid), rownames(G0)); stopifnot(!anyNA(rows))
planted <- rbind(data.frame(rsid = causal$rsid, target = "Y",   eff = causal$b_Y),
                 data.frame(rsid = causal$rsid[causal$cat == "causal_YS"], target = "plt",
                            eff = causal$g_Z[causal$cat == "causal_YS"]),
                 data.frame(rsid = snull$rsid, target = snull$col, eff = snull$g_S))
targets <- c("Y", names(S_COLS))
score <- matrix(0, n, length(targets), dimnames = list(NULL, targets))
usnp <- unique(planted$rsid)
cnt_freq <- setNames(numeric(length(usnp)), usnp)
for (ch in split(usnp, ceiling(seq_along(usnp) / 200))) {
  Gm <- as.matrix(G0[rows, ch, drop = FALSE]); Gm[is.na(Gm)] <- 0
  mu <- colMeans(Gm); cnt_freq[ch] <- mu / 2
  Gm <- sweep(Gm, 2, mu)
  for (tg in targets) {
    pp <- planted[planted$target == tg & planted$rsid %in% ch, ]
    if (nrow(pp)) score[, tg] <- score[, tg] + as.numeric(Gm[, pp$rsid, drop = FALSE] %*% pp$eff)
  }
}
tick("planted scores computed; var: ", paste(sprintf("%s=%.4g", targets, apply(score, 2, var)), collapse = " "))
cat("share of Var(S0) added per column:\n"); print(round(apply(score[, names(S_COLS)], 2, var) / sdS0^2, 4))

snps <- rbind(
  data.frame(causal[, c("rsid", "chr", "pos", "maf", "alt_freq", "call", "cat")], col = NA_character_,
             q = causal$q, qS = ifelse(causal$cat == "causal_YS", QZ_YS, 0), b_Y = causal$b_Y, g_S = causal$g_Z),
  data.frame(snull[, c("rsid", "chr", "pos", "maf", "alt_freq", "call", "cat")], col = snull$col,
             q = 0, qS = snull$qS, b_Y = 0, g_S = snull$g_S),
  data.frame(rnull[, c("rsid", "chr", "pos", "maf", "alt_freq", "call", "cat")], col = NA_character_,
             q = 0, qS = 0, b_Y = 0, g_S = 0))
snps$col[snps$cat == "causal_YS"] <- "plt"
snps$counted_freq <- unname(cnt_freq[snps$rsid])
## b_Y is per BEDMatrix-counted allele (the dosage the scan tests)
saveRDS(list(snps = snps, H_COLS = H_COLS, Z_COL = Z_COL, sdS0 = sdS0, QZ_YS = QZ_YS,
             seeds = c(perm = SEED_PERM, snp = SEED_SNP, eff = SEED_EFF),
             n = n, moved = moved, created = Sys.time(), code = "prep_semisyn.R"),
        paste0(OUT, "design.rds"))
saveRDS(list(perm = perm, score = score, eid = FD$analysis_eid),
        paste0(OUT, "design_indiv.rds"))
file.remove(kf)
tick("saved design.rds (", nrow(snps), " SNPs) and design_indiv.rds")

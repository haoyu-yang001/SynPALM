## =====================================================================
## make_folds_once.R
## Build the five-fold, relatedness-cluster-based assignment ONCE and cache
## it. Every protein job then reads this file, so all proteins share the
## same folds, the same analysis sample, the same GRM submatrix and the same
## set of independent individuals.
##
## Run this once, before launching any protein job.
## =====================================================================

library(data.table)
library(Matrix)

setwd("<PROJECT_DIR>")
source("code/funcs_v2.R")

K          <- 5
FOLD_SEED  <- 20260813        # report this in the Methods
OUT        <- "data/crossfit_folds_K5_afr.rds"

## ---- GRM ------------------------------------------------------------
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData"))
sgrm <- sgrm * 2
colnames(sgrm) <- sapply(colnames(sgrm), function(x) strsplit(x, "_")[[1]][1])
rownames(sgrm) <- sapply(rownames(sgrm), function(x) strsplit(x, "_")[[1]][1])
eids_sGRM <- colnames(sgrm)

## ---- analysis cohort -------------------------------------------------
## Under cross-fitting the analysis set is the SAME for every protein:
## the full post-QC afr cohort. Nothing is held out for training.
df_all_variables_afr_orignal <- fread(
  "data/extracted_pheno_data528_afr.tsv", sep = "\t", header = TRUE, data.table = FALSE)

common_ids   <- eids_sGRM[eids_sGRM %in% df_all_variables_afr_orignal$eid]
analysis_eid <- common_ids                       # this ordering is used everywhere
index_all    <- match(as.character(analysis_eid), eids_sGRM)
GRM          <- sgrm[index_all, index_all]
n_all        <- length(analysis_eid)
cat("analysis cohort n =", n_all, "\n")

## ---- relatedness clusters -------------------------------------------
## find_blocks_vectorized() returns a list of index vectors: the connected
## components of the sparse GRM. These are exactly the relatedness clusters
## we need, and you already use them to pick independent individuals.
blocks <- find_blocks_vectorized(GRM)
sz     <- lengths(blocks)
cat("clusters:", length(blocks),
    " singletons:", sum(sz == 1),
    " largest:", max(sz), "\n")

## ---- assign WHOLE clusters to folds ---------------------------------
## Multi-member clusters first, largest-first into the emptiest fold, then
## the singletons dealt out to equalise the fold sizes. Assigning clusters
## rather than individuals is what makes Yhat_i independent of the
## phenotypes of ALL relatives of i, not just of i.
set.seed(FOLD_SEED)
fold_of_block <- integer(length(blocks))
load_k        <- numeric(K)

multi <- which(sz > 1)
for (b in multi[order(sz[multi], decreasing = TRUE)]) {
  k <- which.min(load_k)
  fold_of_block[b] <- k
  load_k[k] <- load_k[k] + sz[b]
}

sing <- which(sz == 1)
n_s  <- length(sing)
perm <- sample(sing)
target <- sum(sz) / K
need   <- as.integer(floor(pmax(0, target - load_k)))
if (sum(need) > n_s) need <- as.integer(round(need * n_s / sum(need)))
need[K] <- n_s - sum(need[-K])
if (any(need < 0)) {                       # pathological case, fall back
  fold_of_block[perm] <- rep_len(seq_len(K), n_s)
} else {
  fold_of_block[perm] <- rep(seq_len(K), times = need)
}

fold_vec <- integer(n_all)
for (b in seq_along(blocks)) fold_vec[blocks[[b]]] <- fold_of_block[b]
stopifnot(all(fold_vec %in% seq_len(K)))

## ---- sanity checks ---------------------------------------------------
cat("fold sizes:\n"); print(table(fold_vec))
## no relative pairs straddle a fold boundary:
bad <- sum(sapply(blocks, function(ix) length(unique(fold_vec[ix])))> 1)
cat("clusters split across folds (must be 0):", bad, "\n")
stopifnot(bad == 0)

## ---- independent individuals (now protein-invariant) -----------------
independent_indices <- with_seed(25324, {
  sapply(blocks, function(x) x[sample.int(length(x), 1)])
})
cat("independent individuals:", length(independent_indices), "\n")

saveRDS(list(analysis_eid        = analysis_eid,
             index_all           = index_all,
             fold                = fold_vec,
             blocks              = blocks,
             independent_indices = independent_indices,
             K                   = K,
             seed                = FOLD_SEED),
        file = OUT)
cat("written:", OUT, "\n")

## Fold assignments with other random seeds (revision 2026-10-19; reply to Reviewer 1,
## Comment 4): the algorithm of code/synpalm_main_code/make_folds_once.R (whole
## relatedness clusters; multi-member clusters largest-first into the emptiest fold,
## singletons shuffled with the seed and dealt out to equalise the fold sizes),
## applied to the clusters and cohort of data/crossfit_folds_K5.rds (seed 20260813).
##   Rscript code/post_cf3/sens/make_folds_seed.R <seed>
## Writes data/crossfit_folds_K5_seed<seed>.rds (same fields as the main file)
setwd("<PROJECT_DIR>")
FOLD_SEED <- as.integer(commandArgs(TRUE)[1])
FD <- readRDS("data/crossfit_folds_K5.rds")
blocks <- FD$blocks; K <- FD$K; sz <- lengths(blocks); n_all <- length(FD$fold)
set.seed(FOLD_SEED)
fold_of_block <- integer(length(blocks)); load_k <- numeric(K)
multi <- which(sz > 1)
for (b in multi[order(sz[multi], decreasing = TRUE)]) {
  k <- which.min(load_k); fold_of_block[b] <- k; load_k[k] <- load_k[k] + sz[b]
}
sing <- which(sz == 1); n_s <- length(sing); perm <- sample(sing)
target <- sum(sz) / K
need <- as.integer(floor(pmax(0, target - load_k)))
if (sum(need) > n_s) need <- as.integer(round(need * n_s / sum(need)))
need[K] <- n_s - sum(need[-K])
fold_of_block[perm] <- rep(seq_len(K), times = need)
fold_vec <- integer(n_all)
for (b in seq_along(blocks)) fold_vec[blocks[[b]]] <- fold_of_block[b]
stopifnot(all(fold_vec %in% seq_len(K)), sum(sapply(blocks, function(ix) length(unique(fold_vec[ix]))) > 1) == 0)
## check: the main seed reproduces the main assignment
cat("fold sizes:", table(fold_vec), "; agreement with the main folds:", round(mean(fold_vec == FD$fold), 3), "\n")
FD$fold <- fold_vec; FD$seed <- FOLD_SEED
saveRDS(FD, sprintf("data/crossfit_folds_K5_seed%d.rds", FOLD_SEED))

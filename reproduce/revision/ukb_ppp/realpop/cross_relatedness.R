## Relatedness within and across the labelled / unlabelled samples of the 24 proteins
## (revision 2026-10-20; reply to Reviewer 1, Comment 1): pairs of the analysis sample
## with sparse-GRM entry (2 x kinship) >= 0.0884 (kinship >= 0.0442, third degree or
## closer), by labelled-labelled, labelled-unlabelled, unlabelled-unlabelled, and the
## number of labelled individuals with at least one unlabelled relative.
##   Rscript code/post_cf3/realpop/cross_relatedness.R
## Writes results/realpop_cf_noPRS_noPC/cross_relatedness.csv
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
setwd("<PROJECT_DIR>")
sgrm <- get(load("<GRM_DIR>/grmout_nRandomSNPs_0.sGRM.RData")) * 2
FD <- readRDS("data/crossfit_folds_K5.rds")
G  <- sgrm[FD$index_all, FD$index_all]; rm(sgrm)
sm <- as.data.table(summary(as(G, "generalMatrix")))[i < j & x >= 0.0884]
pr <- fread("results/sens_cf_noPRS_noPC/proteins.csv")
hdr <- colnames(fread("data/df_allUKB_eur1.tsv", nrows = 0))[49:2968]
rows <- lapply(seq_len(nrow(pr)), function(k) {
  nm <- hdr[pr$index[k]]
  y <- fread("data/df_allUKB_eur1.tsv", select = c("eid", nm))
  lab <- !is.na(y[[nm]][match(FD$analysis_eid, as.character(y$eid))])
  a <- lab[sm$i]; b <- lab[sm$j]; cross <- a != b
  rel_lab <- unique(c(sm$i[cross & a], sm$j[cross & b]))
  data.table(protein = pr$protein[k], n_lab = sum(lab), pairs_lab_lab = sum(a & b), pairs_lab_unlab = sum(cross),
             pairs_unlab_unlab = sum(!a & !b), first_degree_lab_unlab = sum(cross & sm$x >= 0.354),
             lab_with_unlab_relative = length(rel_lab), share_lab_with_unlab_relative = length(rel_lab) / sum(lab))
})
r <- rbindlist(rows); print(r)
fwrite(r, "results/realpop_cf_noPRS_noPC/cross_relatedness.csv")
